
# *------------------------------------------------------------------
# | FILE NAME: 04_runoff_casestudy_gam.R
# | DATE: 04/25/2025
# | CREATED BY: Connor O'Keefe and Jim Stagge
# *----------------------------------------------------------------


###########################################################################
## Set the Paths
###########################################################################
require(here)

### Path for Data and Output
data_path <- file.path(here(), "data")
#output_path <- "/fs/ess/PAS1921"
output_path <- file.path(here(), "output")

### Set up output folders
write_output_path <- output_path
dir.create(write_output_path, recursive=TRUE, showWarnings = FALSE)

### Set up figure folder
write_figures_path <- file.path(output_path, "figures_04_runoff_casestudy_gam")
dir.create(write_figures_path, recursive=TRUE, showWarnings = FALSE)

###########################################################################
###  Load functions
###########################################################################
require(tidyverse)
require(viridis)

### For Dates
require(lubridate)

### USGS
require(dataRetrieval)

### For GAMS
require(mgcv)
require(gratia)

### To create maps
### Just run this once
#install.packages(c("raster", "sf", "stars", "rnaturalearth", "rnaturalearthhires"))
library(sf)
require(raster)
require(stars)
require(scales)
library(rnaturalearth)
library(rnaturalearthhires)

### Read in all custom functions
sapply(list.files(pattern="[.]R$", path="functions/", full.names=TRUE), source)

select <- dplyr::select



###########################################################################
###  Read in GAGES data
###########################################################################
gages_folder <- file.path(data_path, "gages")

basin_id <- read_csv(file.path(gages_folder, "Dataset1_BasinID/BasinID.txt"))

### Cut to only a few surrounding states
mw_basins <- basin_id %>%
	filter(STATE == "OH" | STATE == "IN" | STATE == "KY" | STATE == "MI" | STATE == "WI" | STATE == "WV") %>%
	as.data.frame()
	
head(mw_basins)

###########################################################################
###  Read in file
###########################################################################
parameterCd <- "00060"  # Discharge
startDate <- "1870-10-01"
endDate <- "2023-09-30"

### Create empty data frame
id_list <- mw_basins$STAID

model_df <- data.frame(STAID = id_list)

### Use to help find gages to investigate (using other code/model)
find_df <- inc_extreme_change %>%
	filter(STATE == "OH") %>%
	filter(season == "annual") %>%
	filter(return_period == "10") %>% 
	filter(change_1950_mm > 1.5) %>% 
	filter(LAT_GAGE > "40.5")
	filter(CLASS == "Ref")


########################################
### Sample list of gauges to look at ###
########################################

### Big Sandy Creek at Rockville, WV (03070500)
	### k <- 46
	### Reference
	### 1909 - 2023
	### Generally decreasing max runoff inc extremes
		### Mostly representative of WV ???

### West Fork River at Butcherville, WV (03058500)
	### k <- 30
	### Non-Ref
	### DA = 489.7 sq km
	### 1915 - 2023 
	### Decreasing max runoff inc extremes

### Johns Creek near Meta, KY (03210000)
	### k <- 190 
	### Reference
	### DA = 146.1 sq km
	### 1941 - 2023
	### Generally decreasing max runoff inc extremes
		### 

### Mill Creek near Coshocton OH (03140000)
	### k <- 106
	### Reference
	### 1936 - 2023
	### Generally increasing max runoff inc extremes
	### Mostly representative of OH ???

### Auglaize River near Fort Jennings OH (04186500) - northwest Ohio
	### k <- 920
	### Non-Ref
	### DA = 858.4 sq km
	### 1921 - 2023
	### Increasing max runoff inc extremes

### West Branch Ontonagon River near Bergland, MI (04036000)
	### k <- 616
	### Non-Ref
	### 1942 - 2023
	### DA = 419.1 sq km
	### Slight decrease in max ruonff inc extremes
	
### Mississinewa River at Marion, IN (03326500)
	### k <- 451
	### Non-Ref
	### DA = 1769.7
	### 1923 - 2023
	### Increasing max runoff inc extremes
		### low point around 1970

### Youngs Creek near Edinburgh, IN (03362000)
	### k <- 550
	### Non-Ref
	### DA = 260.4
	### 1942 - 2023
	### Increasing max runoff inc extremes


#######################################################
### Select Gage for case study                      ###
### Auglaize River near Fort Jennings OH (04186500) ###
#######################################################
k <- 920
	
cat(k)
cat("\n")
	
	### Read in flow from USGS
	#gauge_k <- id_list[[k]]
	gauge_k <- model_df$STAID[k]
	flow_k <- readNWISdv(gauge_k, parameterCd, startDate, endDate)
	flow_k <- renameNWISColumns(flow_k)

	flow_k <- flow_k %>%
	  rename(flow = Flow) %>%
	  rename(date = Date) %>%
	  rename(code = Flow_cd) %>%
	  mutate(units = "cfs") %>%
	  select(-agency_cd)
	  

		### Determine drainage area (sq km)
		drain_area_sqkm <- mw_basins$DRAIN_SQKM[k]
		drain_area_sqkm

	### Convert units of flow (ft3/s) to runoff (mm/day), dividing by drainage area
		flow_k <- flow_k %>% 
			rename(flow_cfs = flow) %>%	
			mutate(flow_mm_day = flow_cfs / drain_area_sqkm / (3280.84)^3 * 10^6 * 86400)

		### Separate the flow
		flow_sep <- process_flow_runoff(flow_k)


		###########################################################################
		###  Fit a GEV model for increasing limb
		###########################################################################
		### Cut to only increasing limb
		inc_df <- flow_sep$inc

		### Add a log10 transform
		inc_df <- inc_df %>%
			mutate(diff_log = log10(diff))
		
		### Calculate mean annual runoff
		#annual_inc_df <- inc_df
		
		### Figure out first and last month
		first_year <- min(inc_df$year)
		last_year <- max(inc_df$year)

		### Set the knots
		knot_loc <- list(jdate=seq(1, 365, length.out = 8), year=seq(1900, 2020, by = 20))
		n_knots <- sapply(knot_loc, FUN = length)


		### Trends in all
		inc_tensor <- gam(list(diff_log
			~ te(jdate, year, bs = c("cc", "tp"), k = c(5, 5)),
			~ te(jdate, year, bs = c("cc", "tp"), k = c(5, 5)),
			~ te(jdate, year, bs = c("cc", "tp"), k = c(5, 5))
		  ),
		  data = inc_df,  #### This distrbution is strictly positive
		  family=gevlss(link = list("identity", "identity", "identity")),
		  select = TRUE,
		  optimizer="efs")



		###########################################################################
		###  Plot the parameters
		###########################################################################

		### Figure out first and last year
		first_year <- min(inc_df$year)
		last_year <- max(inc_df$year)
			first_year
			last_year

		### Create blank data frame for months and years
		newdata <- data.frame(date = seq(as.Date(paste0(first_year, "-01-01")), as.Date(paste0(2023, "-12-31")), by = "1 day")) %>%
		#newdata <- data.frame(date = seq(as.Date(paste0(first_year, "-01-01")), as.Date(paste0(last_year, "-12-31")), by = "1 day")) %>%
			mutate(month = month(date), jdate = yday(date), year = year(date)) %>%
			filter(jdate <= 365)

		### Estimate model parameters at these months and years
		est_df <- pred_gev(mod = inc_tensor, newdata = newdata)

		### Calculate max value per year
		max_value <- est_df %>%
			group_by(year) %>%
			summarize(loc_max = max(est_location, na.rm=TRUE), jdate_max = jdate[which.max(est_location)], loc_min = min(est_location, na.rm=TRUE), jdate_min = jdate[which.min(est_location)])

			### Make plot
			p <- ggplot(est_df %>% filter(year>=first_year & year<= last_year), aes(x=jdate, y=year)) %>%
				+ geom_raster(aes(fill = est_location)) %>%
				#+ geom_path(data = max_value, aes(x=jdate_max), colour = "red") %>%
				#+ geom_path(data = max_value, aes(x=jdate_min), colour = "red", linetype = "dashed") %>%
				+ scale_y_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_x_continuous(name = "Julian Date") %>%
				+ scale_fill_viridis() %>%
				+ theme_bw(14)
			### Save plot
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_est_loc_annual_TENSOR.png")), p,  width = 6.5, height = 4.8, dpi = 600)



			###########################################################################
			###  Find the 2, 5, 10, 25, 50, and 100 year
			###########################################################################
			return_period <- c(2, 5, 10, 25, 50, 100)

			extreme_df <- expand_grid(est_df, data.frame(return_period = return_period)) %>%
				mutate(return_p = 1 - (1/return_period))

			extreme_df <- extreme_df %>%
				rowwise() %>%
				mutate(extreme_inc = evd::qgev(return_p, loc = est_location, scale = est_scale, shape = est_shape)) %>%
			### Reverse the log10 transformed data back to runoff (mm/day)
				mutate(extreme_inc = 10^(extreme_inc))

			### Plot 5-yr daily inc (mm) through calendar year (Jan - Dec)
			p <- ggplot(extreme_df %>% 
					filter(year>=first_year & year<= last_year) %>% 
					filter(year %in% seq(1930,2010, by = 10)) %>% 
					filter(return_period == '5'), 
					aes(x=jdate, y=extreme_inc, color = as.factor(year))) %>%
				+ geom_line() %>%
				+ scale_y_continuous(name = "Daily runofff increase (mm), 5-yr") %>%
				+ scale_x_continuous(name = "Julian Date") %>%
				+ theme_bw(14)
			### Save Plot
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_inc_5yr_calendar_year_line.png")), p,  width = 6.5, height = 4.8, dpi = 600)
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_inc_5yr_calendar_year_line.svg")), p,  width = 6.5, height = 4.8)
			


			### Plot 10-yr daily inc (mm) through calendar year (Jan - Dec)
			p <- ggplot(extreme_df %>% 
					filter(year>=first_year & year<= last_year) %>% 
					filter(year %in% seq(1930,2010, by = 10)) %>% 
					filter(return_period == '10'), 
					aes(x=jdate, y=extreme_inc, color = as.factor(year))) %>%
					#aes(x=jdate, y=extreme_inc)) %>%
				+ geom_line() %>%
				+ scale_y_continuous(name = "Daily runoff (mm/day)") %>%
				+ scale_x_continuous(name = "Julian Date") %>%
				+ theme_bw(14) %>%
				+ theme(legend.position = "none")
			### Save Plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_inc_10yr_calendar_year_line.png")), p,  width = 6.5, height = 4.8, dpi = 1200)
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_inc_10yr_calendar_year_line.png")), p,  width = 6, height = 3, dpi = 1200)
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_inc_10yr_calendar_year_line.svg")), p,  width = 6.5, height = 4.8)



			### Calculate max value per year
			max_value <- extreme_df %>%
				group_by(year, return_period) %>%
				summarize(extreme_max = max(extreme_inc, na.rm=TRUE), jdate_max = jdate[which.max(extreme_inc)]) %>%
				mutate(extreme_max = ifelse(year < first_year, NA, extreme_max)) %>%
				mutate(jdate_max = ifelse(year < first_year, NA, jdate_max)) %>%
				mutate(STAID = model_df$STAID[k])

			### Make plot for each return period
			p <- ggplot(extreme_df %>% filter(year>=1930 & year<=2010), aes(x=jdate, y=year)) %>%
				+ geom_raster(aes(fill = extreme_inc)) %>%
				+ geom_path(data = max_value %>% filter(year>=1930 & year<=2010), aes(x=jdate_max), colour = "red") %>%
				+ scale_y_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_x_continuous(name = "Julian Date") %>%
				+ scale_fill_viridis() %>%
				+ facet_wrap(. ~ return_period) %>%
				+ theme_bw(14)
			### Save plot
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_annual.png")), p,  width = 9, height = 7, dpi = 600)
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_annual.svg")), p,  width = 9, height = 7)



			### Only the 5 year
			p <- ggplot(extreme_df %>% filter(return_period == 5) %>% filter(year>=1930 & year<=2010), aes(x=jdate, y=year)) %>%
				+ geom_raster(aes(fill = extreme_inc)) %>%
				+ geom_path(data = max_value%>% filter(return_period == 5) %>% filter(year>=1930 & year<=2010), aes(x=jdate_max), colour = "red") %>%
				+ scale_y_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_x_continuous(name = "Julian Date") %>%
				+ scale_fill_viridis() %>%
				+ facet_wrap(. ~ return_period) %>%
				+ theme_bw(14)
			### Save plot
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_5yr_annual.png")), p,  width = 6.5, height = 4.8, dpi = 600)
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_5yr_annual.svg")), p,  width = 6.5, height = 4.8)
			
			### Plot annual max single day 5-yr inc through time
			p <- ggplot(max_value %>% filter(return_period == 5) %>% filter(year>=1930 & year<=2010), aes(x=year, y=extreme_max, colour=as.factor(return_period))) %>%
				+ geom_line() %>%
				+ scale_x_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_y_continuous(name = "Max single day runoff increase (mm)") %>%
				+ coord_cartesian(ylim = c(0, 5)) %>%
				#+ coord_cartesian(xlim = c(1920, 2023)) %>%
				+ guides(color = guide_legend(title = "Return Period"))
			### Save plot
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_5yr_annual_line.png")), p,  width = 6.5, height = 4.8, dpi = 600)
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_5yr_annual_line.svg")), p,  width = 6.5, height = 4.8)



			### Only the 10 year
			p <- ggplot(extreme_df %>% filter(return_period == 10) %>% filter(year>=1930 & year<=2010), aes(x=jdate, y=year)) %>%
				+ geom_raster(aes(fill = extreme_inc)) %>%
				+ geom_path(data = max_value%>% filter(return_period == 10) %>% filter(year>=1930 & year<=2010), aes(x=jdate_max), colour = "red") %>%
				+ scale_y_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_x_continuous(name = "Julian Date") %>%
				+ scale_fill_viridis() %>%
				#+ facet_wrap(. ~ return_period) %>%
				+ theme_bw(14) %>%
				+ theme(legend.position = "none")
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_10yr_annual.png")), p,  width = 6.5, height = 4.8, dpi = 1200)
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_10yr_annual.png")), p,  width = 6, height = 6, dpi = 1200)
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_10yr_annual.svg")), p,  width = 6.5, height = 4.8)
			
			### Plot annual max single day 10-yr inc through time
			p <- ggplot(max_value %>% filter(return_period == 10) %>% filter(year>=1930 & year<=2010), aes(x=year, y=extreme_max, colour=as.factor(return_period))) %>%
				+ geom_line() %>%
				+ scale_x_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_y_continuous(name = "Daily Runoff (mm/day)") %>%
				+ coord_cartesian(ylim = c(3, 4.5)) %>% 
				+ theme_bw(14) %>%
				+ theme(legend.position = "none")
				#+ coord_cartesian(xlim = c(1920, 2023)) %>%
				#+ guides(color = guide_legend(title = "Return Period"))
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_10yr_annual_line.png")), p,  width = 6.5, height = 4.8, dpi = 1200)
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_10yr_annual_line.png")), p,  width = 6, height = 3, dpi = 1200)
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_10yr_annual_line.svg")), p,  width = 6.5, height = 4.8)
			
			




			### Plot annual max single day inc extremes (2,5,10,25,50-yr) through time
			p <- ggplot(max_value %>% filter(return_period != 100) %>% filter(year>=1930 & year<=2010), aes(x=year, y=extreme_max, colour=as.factor(return_period))) %>%
				+ geom_line() %>%
				+ scale_x_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_y_continuous(name = "Max single day runoff increase (mm)") %>%
				#+ coord_cartesian(ylim = c(0, 20)) %>%
				+ guides(color = guide_legend(title = "Return Period"))
			### Save plot
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_trends_annual.png")), p,  width = 6.5, height = 4.8, dpi = 600)



			### Calculate percent change in extreme_max relative to 1950
			max_value_perc <- max_value %>%
				group_by(return_period) %>%
				mutate(base_value = extreme_max[year == 1950],
						percent_change = ifelse(is.na(extreme_max) | base_value == 0, NA, ((extreme_max - base_value) / base_value) * 100) ) %>%
						ungroup()

			### Plot change (%) in annual max single day inc extremes (2,5,10,25,50-yr) relative to 1950
			p <- ggplot(max_value_perc %>% filter(return_period != 100) %>% filter(year>=first_year & year<= last_year), aes(x=year, y=percent_change, colour=as.factor(return_period))) %>%
				+ geom_line() %>%
				+ scale_x_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_y_continuous(name = "Change (%) in max single day runoff increase") %>%
				+ coord_cartesian(ylim = c(-100, 100)) %>%
				+ geom_hline(yintercept = 0, linetype = "longdash") %>%
				+ guides(color = guide_legend(title = "Return Period"))
			### Save plot
			ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_trends_perc_annual.png")), p,  width = 6.5, height = 4.8, dpi = 600)


			
			
			

###########################################################################
###  Compare each season individually
###########################################################################
			
			### Filter existing data frame for seasons
			extreme_df_winter <- extreme_df %>%
				filter(month == "1" | month == "2" | month == "12") %>%
				as.data.frame() %>% 
				mutate(year = case_when(month == 12 ~ year+1, 
							TRUE  ~ year))

			extreme_df_spring <- extreme_df %>%
				filter(month == "3" | month == "4" | month == "5") %>%
				as.data.frame()
			
			extreme_df_summer <- extreme_df %>%
				filter(month == "6" | month == "7" | month == "8") %>%
				as.data.frame()

			extreme_df_autumn <- extreme_df %>%
				filter(month == "9" | month == "10" | month == "11") %>%
				as.data.frame()			
				
				
			
			### Filter existing data frame for date in middle of each season
				### Winter: Jan 15
				### Spring: April 15
				### Summer: July 15
				### Autumn: October 15
			extreme_df_winter <- extreme_df %>%
				filter(jdate == "15") %>%
				as.data.frame()
			
			extreme_df_spring <- extreme_df %>%
				filter(jdate == "105") %>%
				as.data.frame()

			extreme_df_summer <- extreme_df %>%
				filter(jdate == "196") %>%
				as.data.frame()				

			extreme_df_autumn <- extreme_df %>%
				filter(jdate == "288") %>%
				as.data.frame()


				
			### Calculate max value per year
			max_value_winter <- extreme_df_winter %>%
				group_by(year, return_period) %>%
				summarize(extreme_max = max(extreme_inc, na.rm=TRUE), jdate_max = jdate[which.max(extreme_inc)]) %>%
				mutate(extreme_max = ifelse(year < first_year, NA, extreme_max)) %>%
				mutate(jdate_max = ifelse(year < first_year, NA, jdate_max)) %>%
				mutate(STAID = model_df$STAID[k])
			
			max_value_spring <- extreme_df_spring %>%
				group_by(year, return_period) %>%
				summarize(extreme_max = max(extreme_inc, na.rm=TRUE), jdate_max = jdate[which.max(extreme_inc)]) %>%
				mutate(extreme_max = ifelse(year < first_year, NA, extreme_max)) %>%
				mutate(jdate_max = ifelse(year < first_year, NA, jdate_max)) %>%
				mutate(STAID = model_df$STAID[k])
				
			max_value_summer <- extreme_df_summer %>%
				group_by(year, return_period) %>%
				summarize(extreme_max = max(extreme_inc, na.rm=TRUE), jdate_max = jdate[which.max(extreme_inc)]) %>%
				mutate(extreme_max = ifelse(year < first_year, NA, extreme_max)) %>%
				mutate(jdate_max = ifelse(year < first_year, NA, jdate_max)) %>%
				mutate(STAID = model_df$STAID[k])
			
			max_value_autumn <- extreme_df_autumn %>%
				group_by(year, return_period) %>%
				summarize(extreme_max = max(extreme_inc, na.rm=TRUE), jdate_max = jdate[which.max(extreme_inc)]) %>%
				mutate(extreme_max = ifelse(year < first_year, NA, extreme_max)) %>%
				mutate(jdate_max = ifelse(year < first_year, NA, jdate_max)) %>%
				mutate(STAID = model_df$STAID[k])



		### Combine and plot seasonal changes
		max_value_season <- max_value_winter %>%
			bind_rows(max_value_spring) %>%
			bind_rows(max_value_summer) %>%
			bind_rows(max_value_autumn)



		### Plot seasonal max single day inc extremes (2,5,10,25,50-yr) through time
		p <- ggplot(max_value_season %>% filter(return_period != 100) %>% filter(year>=first_year & year<= last_year), aes(x=year, y=extreme_max, colour=as.factor(return_period))) %>%
			+ geom_line() %>%
			+ facet_wrap(~factor(jdate_max, c("15", "105", "196", "288"))) %>%
			+ scale_x_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
			+ scale_y_continuous(name = "Max single day runoff increase (mm)") %>%
			+ guides(color = guide_legend(title = "Return Period"))
		### Save plot
		ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_trends_seasonal.png")), p,  width = 9, height = 7, dpi = 600)
		
		
		
		### Calculate percent change in extreme_max relative to 1970
		max_value_season_perc <- max_value_season %>%
			group_by(return_period, jdate_max) %>%
			mutate(base_value = extreme_max[year == 1970],
					percent_change = ifelse(is.na(extreme_max) | base_value == 0, NA, ((extreme_max - base_value) / base_value) * 100) ) %>%
			ungroup()

		
		### Plot change (%) in autumn max single day inc extremes (2,5,10,25,50-yr) relative to 1970
		p <- ggplot(max_value_season_perc %>% filter(return_period != 100) %>% filter(year>=first_year & year<= last_year), aes(x=year, y=percent_change, colour=as.factor(return_period))) %>%
			+ geom_line() %>%
			+ facet_wrap(~factor(jdate_max, c("15", "105", "196", "288"))) %>%
			+ scale_x_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
			+ scale_y_continuous(name = "Change (%) in max single day runoff increase") %>%
			+ coord_cartesian(ylim = c(-100, 100)) %>%
			+ geom_hline(yintercept = 0, linetype = "longdash") %>%
			+ guides(color = guide_legend(title = "Return Period"))
		### Save plot
		ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_trends_perc_seasonal.png")), p,  width = 9, height = 7, dpi = 600)


















		### Winter
			### Make plot for each return period
			p <- ggplot(extreme_df_winter %>% filter(year>=first_year & year<= last_year), aes(x=jdate, y=year)) %>%
				+ geom_raster(aes(fill = extreme_inc)) %>%
				+ geom_path(data = max_value_winter, aes(x=jdate_max), colour = "red") %>%
				+ scale_y_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_x_continuous(name = "Julian Date") %>%
				+ scale_fill_viridis() %>%
				+ facet_wrap(. ~ return_period) %>%
				+ theme_bw(14)
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_winter.png")), p,  width = 9, height = 7, dpi = 600)



			### Only the 5 year
			p <- ggplot(extreme_df_winter %>% filter(return_period == 5) %>% filter(year>=first_year & year<= last_year), aes(x=jdate, y=year)) %>%
				+ geom_raster(aes(fill = extreme_inc)) %>%
				+ geom_path(data = max_value_winter%>% filter(return_period == 5) %>% filter(year>=first_year & year<= last_year), aes(x=jdate_max), colour = "red") %>%
				+ scale_y_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_x_continuous(name = "Julian Date") %>%
				+ scale_fill_viridis() %>%
				+ facet_wrap(. ~ return_period) %>%
				+ theme_bw(14)
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_5yr_winter.png")), p,  width = 6.5, height = 4.8, dpi = 600)



			### Plot winter max single day inc extremes (2,5,10,25,50-yr) through time
			p <- ggplot(max_value_winter %>% filter(return_period != 100) %>% filter(year>=first_year & year<= last_year), aes(x=year, y=extreme_max, colour=as.factor(return_period))) %>%
				+ geom_line() %>%
				+ scale_x_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_y_continuous(name = "Max single day runoff increase (mm)") %>%
				+ guides(color = guide_legend(title = "Return Period"))
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_trends_winter.png")), p,  width = 6.5, height = 4.8, dpi = 600)



			### Calculate percent change in extreme_max relative to 1950
			max_value_winter_perc <- max_value_winter %>%
				group_by(return_period) %>%
				mutate(base_value = extreme_max[year == 1970],
						percent_change = ifelse(is.na(extreme_max) | base_value == 0, NA, ((extreme_max - base_value) / base_value) * 100) ) %>%
						ungroup()

			### Plot change (%) in winter max single day inc extremes (2,5,10,25,50-yr) relative to 1970
			p <- ggplot(max_value_winter_perc %>% filter(return_period != 100) %>% filter(year>=first_year & year<= last_year), aes(x=year, y=percent_change, colour=as.factor(return_period))) %>%
				+ geom_line() %>%
				+ scale_x_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_y_continuous(name = "Change (%) in max single day runoff increase") %>%
				+ coord_cartesian(ylim = c(-100, 100)) %>%
				+ geom_hline(yintercept = 0, linetype = "longdash") %>%
				+ guides(color = guide_legend(title = "Return Period"))
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_trends_perc_winter.png")), p,  width = 6.5, height = 4.8, dpi = 600)




		### Spring
			### Make plot for each return period
			p <- ggplot(extreme_df_spring %>% filter(year>=first_year & year<= last_year), aes(x=jdate, y=year)) %>%
				+ geom_raster(aes(fill = extreme_inc)) %>%
				+ geom_path(data = max_value_spring, aes(x=jdate_max), colour = "red") %>%
				+ scale_y_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_x_continuous(name = "Julian Date") %>%
				+ scale_fill_viridis() %>%
				+ facet_wrap(. ~ return_period) %>%
				+ theme_bw(14)
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_spring.png")), p,  width = 9, height = 7, dpi = 600)



			### Only the 5 year
			p <- ggplot(extreme_df_spring %>% filter(return_period == 5) %>% filter(year>=first_year & year<= last_year), aes(x=jdate, y=year)) %>%
				+ geom_raster(aes(fill = extreme_inc)) %>%
				+ geom_path(data = max_value_spring%>% filter(return_period == 5) %>% filter(year>=first_year & year<= last_year), aes(x=jdate_max), colour = "red") %>%
				+ scale_y_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_x_continuous(name = "Julian Date") %>%
				+ scale_fill_viridis() %>%
				+ facet_wrap(. ~ return_period) %>%
				+ theme_bw(14)
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_5yr_spring.png")), p,  width = 6.5, height = 4.8, dpi = 600)



			### Plot spring max single day inc extremes (2,5,10,25,50-yr) through time
			p <- ggplot(max_value_spring %>% filter(return_period != 100) %>% filter(year>=first_year & year<= last_year), aes(x=year, y=extreme_max, colour=as.factor(return_period))) %>%
				+ geom_line() %>%
				+ scale_x_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_y_continuous(name = "Max single day runoff increase (mm)") %>%
				+ guides(color = guide_legend(title = "Return Period"))
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_trends_spring.png")), p,  width = 6.5, height = 4.8, dpi = 600)



			### Calculate percent change in extreme_max relative to 1950
			max_value_spring_perc <- max_value_spring %>%
				group_by(return_period) %>%
				mutate(base_value = extreme_max[year == 1970],
						percent_change = ifelse(is.na(extreme_max) | base_value == 0, NA, ((extreme_max - base_value) / base_value) * 100) ) %>%
						ungroup()

			### Plot change (%) in spring max single day inc extremes (2,5,10,25,50-yr) relative to 1970
			p <- ggplot(max_value_spring_perc %>% filter(return_period != 100) %>% filter(year>=first_year & year<= last_year), aes(x=year, y=percent_change, colour=as.factor(return_period))) %>%
				+ geom_line() %>%
				+ scale_x_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_y_continuous(name = "Change (%) in max single day runoff increase") %>%
				+ coord_cartesian(ylim = c(-100, 100)) %>%
				+ geom_hline(yintercept = 0, linetype = "longdash") %>%
				+ guides(color = guide_legend(title = "Return Period"))
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_trends_perc_spring.png")), p,  width = 6.5, height = 4.8, dpi = 600)
	




		### Summer
			### Make plot for each return period
			p <- ggplot(extreme_df_summer %>% filter(year>=first_year & year<= last_year), aes(x=jdate, y=year)) %>%
				+ geom_raster(aes(fill = extreme_inc)) %>%
				+ geom_path(data = max_value_summer, aes(x=jdate_max), colour = "red") %>%
				+ scale_y_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_x_continuous(name = "Julian Date") %>%
				+ scale_fill_viridis() %>%
				+ facet_wrap(. ~ return_period) %>%
				+ theme_bw(14)
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_summer.png")), p,  width = 9, height = 7, dpi = 600)



			### Only the 5 year
			p <- ggplot(extreme_df_summer %>% filter(return_period == 5) %>% filter(year>=first_year & year<= last_year), aes(x=jdate, y=year)) %>%
				+ geom_raster(aes(fill = extreme_inc)) %>%
				+ geom_path(data = max_value_summer%>% filter(return_period == 5) %>% filter(year>=first_year & year<= last_year), aes(x=jdate_max), colour = "red") %>%
				+ scale_y_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_x_continuous(name = "Julian Date") %>%
				+ scale_fill_viridis() %>%
				+ facet_wrap(. ~ return_period) %>%
				+ theme_bw(14)
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_5yr_summer.png")), p,  width = 6.5, height = 4.8, dpi = 600)



			### Plot summer max single day inc extremes (2,5,10,25,50-yr) through time
			p <- ggplot(max_value_summer %>% filter(return_period != 100) %>% filter(year>=first_year & year<= last_year), aes(x=year, y=extreme_max, colour=as.factor(return_period))) %>%
				+ geom_line() %>%
				+ scale_x_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_y_continuous(name = "Max single day runoff increase (mm)") %>%
				+ guides(color = guide_legend(title = "Return Period"))
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_trends_summer.png")), p,  width = 6.5, height = 4.8, dpi = 600)



			### Calculate percent change in extreme_max relative to 1950
			max_value_summer_perc <- max_value_summer %>%
				group_by(return_period) %>%
				mutate(base_value = extreme_max[year == 1970],
						percent_change = ifelse(is.na(extreme_max) | base_value == 0, NA, ((extreme_max - base_value) / base_value) * 100) ) %>%
						ungroup()

			### Plot change (%) in summer max single day inc extremes (2,5,10,25,50-yr) relative to 1970
			p <- ggplot(max_value_summer_perc %>% filter(return_period != 100) %>% filter(year>=first_year & year<= last_year), aes(x=year, y=percent_change, colour=as.factor(return_period))) %>%
				+ geom_line() %>%
				+ scale_x_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_y_continuous(name = "Change (%) in max single day runoff increase") %>%
				+ coord_cartesian(ylim = c(-100, 100)) %>%
				+ geom_hline(yintercept = 0, linetype = "longdash") %>%
				+ guides(color = guide_legend(title = "Return Period"))
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_trends_perc_summer.png")), p,  width = 6.5, height = 4.8, dpi = 600)





		### Autumn
			### Make plot for each return period
			p <- ggplot(extreme_df_autumn %>% filter(year>=first_year & year<= last_year), aes(x=jdate, y=year)) %>%
				+ geom_raster(aes(fill = extreme_inc)) %>%
				+ geom_path(data = max_value_autumn, aes(x=jdate_max), colour = "red") %>%
				+ scale_y_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_x_continuous(name = "Julian Date") %>%
				+ scale_fill_viridis() %>%
				+ facet_wrap(. ~ return_period) %>%
				+ theme_bw(14)
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_autumn.png")), p,  width = 9, height = 7, dpi = 600)



			### Only the 5 year
			p <- ggplot(extreme_df_autumn %>% filter(return_period == 5) %>% filter(year>=first_year & year<= last_year), aes(x=jdate, y=year)) %>%
				+ geom_raster(aes(fill = extreme_inc)) %>%
				+ geom_path(data = max_value_autumn%>% filter(return_period == 5) %>% filter(year>=first_year & year<= last_year), aes(x=jdate_max), colour = "red") %>%
				+ scale_y_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_x_continuous(name = "Julian Date") %>%
				+ scale_fill_viridis() %>%
				+ facet_wrap(. ~ return_period) %>%
				+ theme_bw(14)
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_5yr_autumn.png")), p,  width = 6.5, height = 4.8, dpi = 600)



			### Plot autumn max single day inc extremes (2,5,10,25,50-yr) through time
			p <- ggplot(max_value_autumn %>% filter(return_period != 100) %>% filter(year>=first_year & year<= last_year), aes(x=year, y=extreme_max, colour=as.factor(return_period))) %>%
				+ geom_line() %>%
				+ scale_x_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_y_continuous(name = "Max single day runoff increase (mm)") %>%
				+ guides(color = guide_legend(title = "Return Period"))
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_trends_autumn.png")), p,  width = 6.5, height = 4.8, dpi = 600)



			### Calculate percent change in extreme_max relative to 1950
			max_value_autumn_perc <- max_value_autumn %>%
				group_by(return_period) %>%
				mutate(base_value = extreme_max[year == 1970],
						percent_change = ifelse(is.na(extreme_max) | base_value == 0, NA, ((extreme_max - base_value) / base_value) * 100) ) %>%
						ungroup()

			### Plot change (%) in autumn max single day inc extremes (2,5,10,25,50-yr) relative to 1970
			p <- ggplot(max_value_autumn_perc %>% filter(return_period != 100) %>% filter(year>=first_year & year<= last_year), aes(x=year, y=percent_change, colour=as.factor(return_period))) %>%
				+ geom_line() %>%
				+ scale_x_continuous(name = "Year", breaks = seq(1800,2040, by = 10)) %>%
				+ scale_y_continuous(name = "Change (%) in max single day runoff increase") %>%
				+ coord_cartesian(ylim = c(-100, 100)) %>%
				+ geom_hline(yintercept = 0, linetype = "longdash") %>%
				+ guides(color = guide_legend(title = "Return Period"))
			### Save plot
			#ggsave(file.path(write_figures_path, paste0(gauge_k, "_extreme_inc_trends_perc_autumn.png")), p,  width = 6.5, height = 4.8, dpi = 600)
	
			
			
			
