
# *------------------------------------------------------------------
# | FILE NAME: 03_runoff_daily_gam.R
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
write_figures_path <- file.path(output_path, "figures_03_runoff_daily_gam")
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
###  Read in Gauge data
###########################################################################
mw_basins <- read_csv(paste0(data_path, '/seasonal_flow_trends.csv')) 

mw_basins <- mw_basins %>%
	select(STAID, STANAME, DRAIN_SQKM, HUC02, LAT_GAGE, LNG_GAGE, STATE, CLASS, AGGECOREGION)

head(mw_basins)



###  Read in flow
parameterCd <- "00060"  # Discharge
startDate <- "1870-10-01"
endDate <- "2023-09-30"

### Create empty data frame
id_list <- mw_basins$STAID

model_df <- data.frame(STAID = id_list)


### Run loop for each gauge
for (k in seq(1, dim(model_df)[1])){

cat(k)
cat("\n")
	
	### Read in flow from USGS
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

	### Figure out first and last month
	first_year <- min(inc_df$year)
	last_year <- max(inc_df$year)

	### Create blank data frame for months and years
	newdata <- data.frame(date = seq(as.Date(paste0(1920, "-01-01")), as.Date(paste0(2023, "-12-31")), by = "1 day")) %>%
		mutate(month = month(date), jdate = yday(date), year = year(date)) %>%
		filter(jdate <= 365)

	### Estimate model parameters at these months and years
	est_df <- pred_gev(mod = inc_tensor, newdata = newdata)



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



		### Calculate max value per year
		max_value <- extreme_df %>%
			group_by(year, return_period) %>%
			summarize(extreme_max = max(extreme_inc, na.rm=TRUE), jdate_max = jdate[which.max(extreme_inc)]) %>%
			mutate(extreme_max = ifelse(year < first_year, NA, extreme_max)) %>%
			mutate(jdate_max = ifelse(year < first_year, NA, jdate_max)) %>%
			mutate(STAID = model_df$STAID[k])
			
		max_value_df_temp <- max_value



		### Calculate change since 1920 (1920 - 2010)
		change_1920_2yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 2]) - (max_value$extreme_max[max_value$year == 1920 & max_value$return_period == 2])
		change_1920_2yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 2]) - (max_value$jdate_max[max_value$year == 1920 & max_value$return_period == 2])
		
		change_1920_5yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 5]) - (max_value$extreme_max[max_value$year == 1920 & max_value$return_period == 5])
		change_1920_5yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 5]) - (max_value$jdate_max[max_value$year == 1920 & max_value$return_period == 5])
		
		change_1920_10yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 10]) - (max_value$extreme_max[max_value$year == 1920 & max_value$return_period == 10])
		change_1920_10yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 10]) - (max_value$jdate_max[max_value$year == 1920 & max_value$return_period == 10])
		
		change_1920_25yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 25]) - (max_value$extreme_max[max_value$year == 1920 & max_value$return_period == 25])
		change_1920_25yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 25]) - (max_value$jdate_max[max_value$year == 1920 & max_value$return_period == 25])
		
		change_1920_50yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 50]) - (max_value$extreme_max[max_value$year == 1920 & max_value$return_period == 50])
		change_1920_50yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 50]) - (max_value$jdate_max[max_value$year == 1920 & max_value$return_period == 50])
		
		change_1920_100yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 100]) - (max_value$extreme_max[max_value$year == 1920 & max_value$return_period == 100])
		change_1920_100yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 100]) - (max_value$jdate_max[max_value$year == 1920 & max_value$return_period == 100])



		### Calculate change since 1950 (1950 - 2010)
		change_1950_2yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 2]) - (max_value$extreme_max[max_value$year == 1950 & max_value$return_period == 2])
		change_1950_2yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 2]) - (max_value$jdate_max[max_value$year == 1950 & max_value$return_period == 2])
		
		change_1950_5yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 5]) - (max_value$extreme_max[max_value$year == 1950 & max_value$return_period == 5])
		change_1950_5yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 5]) - (max_value$jdate_max[max_value$year == 1950 & max_value$return_period == 5])
		
		change_1950_10yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 10]) - (max_value$extreme_max[max_value$year == 1950 & max_value$return_period == 10])
		change_1950_10yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 10]) - (max_value$jdate_max[max_value$year == 1950 & max_value$return_period == 10])
		
		change_1950_25yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 25]) - (max_value$extreme_max[max_value$year == 1950 & max_value$return_period == 25])
		change_1950_25yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 25]) - (max_value$jdate_max[max_value$year == 1950 & max_value$return_period == 25])
		
		change_1950_50yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 50]) - (max_value$extreme_max[max_value$year == 1950 & max_value$return_period == 50])
		change_1950_50yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 50]) - (max_value$jdate_max[max_value$year == 1950 & max_value$return_period == 50])
		
		change_1950_100yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 100]) - (max_value$extreme_max[max_value$year == 1950 & max_value$return_period == 100])
		change_1950_100yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 100]) - (max_value$jdate_max[max_value$year == 1950 & max_value$return_period == 100])
	
		

		### Calculate change since 1970 (1970 - 2010)
		change_1970_2yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 2]) - (max_value$extreme_max[max_value$year == 1970 & max_value$return_period == 2])
		change_1970_2yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 2]) - (max_value$jdate_max[max_value$year == 1970 & max_value$return_period == 2])
		
		change_1970_5yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 5]) - (max_value$extreme_max[max_value$year == 1970 & max_value$return_period == 5])
		change_1970_5yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 5]) - (max_value$jdate_max[max_value$year == 1970 & max_value$return_period == 5])
		
		change_1970_10yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 10]) - (max_value$extreme_max[max_value$year == 1970 & max_value$return_period == 10])
		change_1970_10yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 10]) - (max_value$jdate_max[max_value$year == 1970 & max_value$return_period == 10])
		
		change_1970_25yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 25]) - (max_value$extreme_max[max_value$year == 1970 & max_value$return_period == 25])
		change_1970_25yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 25]) - (max_value$jdate_max[max_value$year == 1970 & max_value$return_period == 25])
		
		change_1970_50yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 50]) - (max_value$extreme_max[max_value$year == 1970 & max_value$return_period == 50])
		change_1970_50yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 50]) - (max_value$jdate_max[max_value$year == 1970 & max_value$return_period == 50])
		
		change_1970_100yr <- (max_value$extreme_max[max_value$year == 2010 & max_value$return_period == 100]) - (max_value$extreme_max[max_value$year == 1970 & max_value$return_period == 100])
		change_1970_100yr_date <- (max_value$jdate_max[max_value$year == 2010 & max_value$return_period == 100]) - (max_value$jdate_max[max_value$year == 1970 & max_value$return_period == 100])


		
		### Create empty dataframe to store calculated changes
		max_value_change_annual <- data.frame(return_period = return_period) %>%
			mutate(change_1920_mm = NA) %>%
			mutate(change_1920_mm_perc = NA) %>%
			mutate(change_1920_date = NA) %>%
			mutate(change_1950_mm = NA) %>%
			mutate(change_1950_mm_perc = NA) %>%
			mutate(change_1950_date = NA) %>%
			mutate(change_1970_mm = NA) %>%
			mutate(change_1970_mm_perc = NA) %>%
			mutate(change_1970_date = NA) %>%
			mutate(season = 'annual') %>%
			mutate(STAID = model_df$STAID[k]) 

			### Store values
			max_value_change_annual$change_1920_mm[return_period == 2] <- change_1920_2yr
			max_value_change_annual$change_1920_mm[return_period == 5] <- change_1920_5yr
			max_value_change_annual$change_1920_mm[return_period == 10] <- change_1920_10yr
			max_value_change_annual$change_1920_mm[return_period == 25] <- change_1920_25yr
			max_value_change_annual$change_1920_mm[return_period == 50] <- change_1920_50yr
			max_value_change_annual$change_1920_mm[return_period == 100] <- change_1920_100yr

			max_value_change_annual$change_1920_mm_perc[return_period == 2] <- (change_1920_2yr / (max_value$extreme_max[max_value$year == 1920 & max_value$return_period == 2])) * 100
			max_value_change_annual$change_1920_mm_perc[return_period == 5] <- (change_1920_5yr / (max_value$extreme_max[max_value$year == 1920 & max_value$return_period == 5])) * 100
			max_value_change_annual$change_1920_mm_perc[return_period == 10] <- (change_1920_10yr / (max_value$extreme_max[max_value$year == 1920 & max_value$return_period == 10])) * 100
			max_value_change_annual$change_1920_mm_perc[return_period == 25] <- (change_1920_25yr / (max_value$extreme_max[max_value$year == 1920 & max_value$return_period == 25])) * 100
			max_value_change_annual$change_1920_mm_perc[return_period == 50] <- (change_1920_50yr / (max_value$extreme_max[max_value$year == 1920 & max_value$return_period == 50])) * 100
			max_value_change_annual$change_1920_mm_perc[return_period == 100] <- (change_1920_100yr / (max_value$extreme_max[max_value$year == 1920 & max_value$return_period == 100])) * 100

			
			max_value_change_annual$change_1920_date[return_period == 2] <- change_1920_2yr_date
			max_value_change_annual$change_1920_date[return_period == 5] <- change_1920_5yr_date
			max_value_change_annual$change_1920_date[return_period == 10] <- change_1920_10yr_date
			max_value_change_annual$change_1920_date[return_period == 25] <- change_1920_25yr_date
			max_value_change_annual$change_1920_date[return_period == 50] <- change_1920_50yr_date
			max_value_change_annual$change_1920_date[return_period == 100] <- change_1920_100yr_date
							
			max_value_change_annual$change_1950_mm[return_period == 2] <- change_1950_2yr
			max_value_change_annual$change_1950_mm[return_period == 5] <- change_1950_5yr
			max_value_change_annual$change_1950_mm[return_period == 10] <- change_1950_10yr
			max_value_change_annual$change_1950_mm[return_period == 25] <- change_1950_25yr
			max_value_change_annual$change_1950_mm[return_period == 50] <- change_1950_50yr
			max_value_change_annual$change_1950_mm[return_period == 100] <- change_1950_100yr

			max_value_change_annual$change_1950_mm_perc[return_period == 2] <- (change_1950_2yr / (max_value$extreme_max[max_value$year == 1950 & max_value$return_period == 2])) * 100
			max_value_change_annual$change_1950_mm_perc[return_period == 5] <- (change_1950_5yr / (max_value$extreme_max[max_value$year == 1950 & max_value$return_period == 5])) * 100
			max_value_change_annual$change_1950_mm_perc[return_period == 10] <- (change_1950_10yr / (max_value$extreme_max[max_value$year == 1950 & max_value$return_period == 10])) * 100
			max_value_change_annual$change_1950_mm_perc[return_period == 25] <- (change_1950_25yr / (max_value$extreme_max[max_value$year == 1950 & max_value$return_period == 25])) * 100
			max_value_change_annual$change_1950_mm_perc[return_period == 50] <- (change_1950_50yr / (max_value$extreme_max[max_value$year == 1950 & max_value$return_period == 50])) * 100
			max_value_change_annual$change_1950_mm_perc[return_period == 100] <- (change_1950_100yr / (max_value$extreme_max[max_value$year == 1950 & max_value$return_period == 100])) * 100
			
			max_value_change_annual$change_1950_date[return_period == 2] <- change_1950_2yr_date
			max_value_change_annual$change_1950_date[return_period == 5] <- change_1950_5yr_date
			max_value_change_annual$change_1950_date[return_period == 10] <- change_1950_10yr_date
			max_value_change_annual$change_1950_date[return_period == 25] <- change_1950_25yr_date
			max_value_change_annual$change_1950_date[return_period == 50] <- change_1950_50yr_date
			max_value_change_annual$change_1950_date[return_period == 100] <- change_1950_100yr_date
				
			max_value_change_annual$change_1970_mm[return_period == 2] <- change_1970_2yr
			max_value_change_annual$change_1970_mm[return_period == 5] <- change_1970_5yr
			max_value_change_annual$change_1970_mm[return_period == 10] <- change_1970_10yr
			max_value_change_annual$change_1970_mm[return_period == 25] <- change_1970_25yr
			max_value_change_annual$change_1970_mm[return_period == 50] <- change_1970_50yr
			max_value_change_annual$change_1970_mm[return_period == 100] <- change_1970_100yr

			max_value_change_annual$change_1970_mm_perc[return_period == 2] <- (change_1970_2yr / (max_value$extreme_max[max_value$year == 1970 & max_value$return_period == 2])) * 100
			max_value_change_annual$change_1970_mm_perc[return_period == 5] <- (change_1970_5yr / (max_value$extreme_max[max_value$year == 1970 & max_value$return_period == 5])) * 100
			max_value_change_annual$change_1970_mm_perc[return_period == 10] <- (change_1970_10yr / (max_value$extreme_max[max_value$year == 1970 & max_value$return_period == 10])) * 100
			max_value_change_annual$change_1970_mm_perc[return_period == 25] <- (change_1970_25yr / (max_value$extreme_max[max_value$year == 1970 & max_value$return_period == 25])) * 100
			max_value_change_annual$change_1970_mm_perc[return_period == 50] <- (change_1970_50yr / (max_value$extreme_max[max_value$year == 1970 & max_value$return_period == 50])) * 100
			max_value_change_annual$change_1970_mm_perc[return_period == 100] <- (change_1970_100yr / (max_value$extreme_max[max_value$year == 1970 & max_value$return_period == 100])) * 100
			
			max_value_change_annual$change_1970_date[return_period == 2] <- change_1970_2yr_date
			max_value_change_annual$change_1970_date[return_period == 5] <- change_1970_5yr_date
			max_value_change_annual$change_1970_date[return_period == 10] <- change_1970_10yr_date
			max_value_change_annual$change_1970_date[return_period == 25] <- change_1970_25yr_date
			max_value_change_annual$change_1970_date[return_period == 50] <- change_1970_50yr_date
			max_value_change_annual$change_1970_date[return_period == 100] <- change_1970_100yr_date

		
		
		### Filter to actual record period (first_year, last_year) and store data frame (could be useful later?)
		max_value <- max_value %>%
		filter(year >= first_year) %>%
		filter(year <= last_year)

		max_value_df_temp <- max_value



		if (k == 1){	#first pass through loop
			max_value_df <- max_value_df_temp
			
			max_value_change <- max_value_change_annual
			
			
			} else {
				max_value_df <- max_value_df %>%
					bind_rows(max_value_df_temp)
					
				max_value_change <- max_value_change %>%
					bind_rows(max_value_change_annual)
				
			}
		
	
	
}





### Combine dataframes with mw_basins
inc_extreme_complete <- mw_basins %>% full_join(max_value_df, by = "STAID")
	inc_extreme_complete <- inc_extreme_complete
	
inc_extreme_change <- mw_basins %>% full_join(max_value_change, by = "STAID")
	inc_extreme_change <- inc_extreme_change

### Write to file
write.csv(inc_extreme_complete, paste0(data_path, '/runoff_extremes_complete.csv')) 
write.csv(inc_extreme_change, paste0(data_path, '/runoff_extremes_change.csv'))

### Create an sf object of gauge locations
new_sf = st_as_sf(inc_extreme_change, coords = c("LNG_GAGE", "LAT_GAGE"), crs = 4326 )
	+ scale_colour_manual(name = "Change", values = c("blue", "black",  "red"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
	+ theme_bw(12)

### Map data comes from https://www.naturalearthdata.com/downloads/
### Download state borders
usa_states <- ne_states(country = 'United States of America', returnclass = 'sf')
ca_states <- ne_states(country = 'Canada', returnclass = 'sf')
mx_states <- ne_states(country = 'Mexico', returnclass = 'sf')
### Download lakes and rivers
lakes <- ne_download(scale = 50, type = 'lakes', category = 'physical', returnclass = 'sf')
rivers <- ne_download(scale = 10, type = 'rivers_lake_centerlines', category = 'physical', returnclass = 'sf')



############################################## ##########################################
### Plot map of changes in single-day inc extremes and DATE (1920, 1950, 1970 - 2010) ###
#########################################################################################

##### 2yr single day inc event: mm/day and date #####
### 1920 - 2010
max_value <- 0.25
max_value_date <- 60

p <- ggplot(data = new_sf %>% drop_na(change_1920_mm) %>% filter(return_period == "2") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1920_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "2yr single-day inc \nchange 1920-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	##+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_2yr_inc_mm_1920_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1920_date) %>% filter(return_period == "2") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1920_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "2yr single-day inc \nchange 1920-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_2yr_inc_date_1920_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)



### 1950 - 2010
p <- ggplot(data = new_sf %>% drop_na(change_1950_mm) %>% filter(return_period == "2") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1950_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "2yr single-day inc \nchange 1950-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_2yr_inc_mm_1950_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1950_date) %>% filter(return_period == "2") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1950_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "2yr single-day inc \nchange 1950-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_2yr_inc_date_1950_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)



### 1970 - 2010
p <- ggplot(data = new_sf %>% drop_na(change_1970_mm) %>% filter(return_period == "2") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1970_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "2yr single-day inc \nchange 1970-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_2yr_inc_mm_1970_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1970_date) %>% filter(return_period == "2") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1970_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "2yr single-day inc \nchange 1970-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_2yr_inc_date_1970_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)



##### 5yr single day inc event: mm/day and date #####
max_value <- 1
max_value_date <- 60

### 1920 - 2010
p <- ggplot(data = new_sf %>% drop_na(change_1920_mm) %>% filter(return_period == "5") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1920_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "5yr single-day inc \nchange 1920-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_5yr_inc_mm_1920_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1920_date) %>% filter(return_period == "5") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1920_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "5yr single-day inc \nchange 1920-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_5yr_inc_date_1920_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)



### 1950 - 2010
p <- ggplot(data = new_sf %>% drop_na(change_1950_mm) %>% filter(return_period == "5") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1950_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "5yr single-day inc \nchange 1950-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_5yr_inc_mm_1950_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1950_date) %>% filter(return_period == "5") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1950_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "5yr single-day inc \nchange 1950-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_5yr_inc_date_1950_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)



### 1970 - 2010
p <- ggplot(data = new_sf %>% drop_na(change_1970_mm) %>% filter(return_period == "5") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1970_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "5yr single-day inc \nchange 1970-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_5yr_inc_mm_1970_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1970_date) %>% filter(return_period == "5") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1970_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "5yr single-day inc \nchange 1970-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_5yr_inc_date_1970_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)



##### 10yr single day inc event: mm/day and date #####
max_value <- 2
max_value_date <- 60

### 1920 - 2010

p <- ggplot(data = new_sf %>% drop_na(change_1920_mm) %>% filter(return_period == "10") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1920_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "10yr single-day inc \nchange 1920-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_10yr_inc_mm_1920_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1920_date) %>% filter(return_period == "10") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1920_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "10yr single-day inc \nchange 1920-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_10yr_inc_date_1920_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)



### 1950 - 2010
p <- ggplot(data = new_sf %>% drop_na(change_1950_mm) %>% filter(return_period == "10") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1950_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "10yr single-day inc \nchange 1950-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_10yr_inc_mm_1950_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1950_date) %>% filter(return_period == "10") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1950_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "10yr single-day inc \nchange 1950-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_10yr_inc_date_1950_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	


### 1970 - 2010
p <- ggplot(data = new_sf %>% drop_na(change_1970_mm) %>% filter(return_period == "10") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1970_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "10yr single-day inc \nchange 1970-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_10yr_inc_mm_1970_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1970_date) %>% filter(return_period == "10") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1970_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "10yr single-day inc \nchange 1970-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_10yr_inc_date_1970_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)



##### 25yr single day inc event: mm/day and date #####
max_value <- 4
max_value_date <- 60

### 1920 - 2010
p <- ggplot(data = new_sf %>% drop_na(change_1920_mm) %>% filter(return_period == "25") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1920_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "25yr single-day inc \nchange 1920-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_25yr_inc_mm_1920_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1920_date) %>% filter(return_period == "25") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1920_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "25yr single-day inc \nchange 1920-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_25yr_inc_date_1920_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)



### 1950 - 2010
p <- ggplot(data = new_sf %>% drop_na(change_1950_mm) %>% filter(return_period == "25") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1950_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "25yr single-day inc \nchange 1950-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_25yr_inc_mm_1950_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1950_date) %>% filter(return_period == "25") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1950_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "25yr single-day inc \nchange 1950-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_25yr_inc_date_1950_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)



### 1970 - 2010
p <- ggplot(data = new_sf %>% drop_na(change_1970_mm) %>% filter(return_period == "25") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1970_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "25yr single-day inc \nchange 1970-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_25yr_inc_mm_1970_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1970_date) %>% filter(return_period == "25") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1970_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "25yr single-day inc \nchange 1970-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_25yr_inc_date_1970_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)





##### 50yr single day inc event: mm/day and date #####
max_value <- 20
max_value_date <- 60

### 1920 - 2010
p <- ggplot(data = new_sf %>% drop_na(change_1920_mm) %>% filter(return_period == "50") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1920_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "50yr single-day inc \nchange 1920-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_50yr_inc_mm_1920_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1920_date) %>% filter(return_period == "50") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1920_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "50yr single-day inc \nchange 1920-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_50yr_inc_date_1920_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)



### 1950 - 2010
p <- ggplot(data = new_sf %>% drop_na(change_1950_mm) %>% filter(return_period == "50") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1950_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "50yr single-day inc \nchange 1950-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_50yr_inc_mm_1950_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1950_date) %>% filter(return_period == "50") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1950_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "50yr single-day inc \nchange 1950-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_50yr_inc_date_1950_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)



### 1970 - 2010
p <- ggplot(data = new_sf %>% drop_na(change_1970_mm) %>% filter(return_period == "50") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1970_mm, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "50yr single-day inc \nchange 1970-2010 \nmm/day", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_50yr_inc_mm_1970_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)

p <- ggplot(data = new_sf %>% drop_na(change_1970_date) %>% filter(return_period == "50") %>% filter(season == "annual")) %>%
	#+ geom_stars(data = background, alpha = 0.7) %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1970_date, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "50yr single-day inc \nchange 1970-2010 \nJulian date", palette = "PuOr", limits = c(-max_value_date, max_value_date), oob=squish, direction = -1) %>%
	+ scale_fill_identity() %>%
	#+ facet_wrap(~factor(season, c("annual"))) %>%
	+ facet_wrap(. ~ season) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_50yr_inc_date_1970_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)





###################################################################################
### Plot map of single-day inc extremes as a % change (1920, 1950, 1970 - 2010) ###
###################################################################################
max_value <- 100

### 1920 - 2010 ###
p <- ggplot(data = new_sf %>% drop_na(change_1920_mm_perc) %>% filter(season == 'annual'))  %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1920_mm_perc, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "Single day inc extremes \n% change 1920-2010 \n Annual", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	+ facet_wrap(. ~ return_period) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_extreme_inc_mm_perc_1920_annual.png"), p,  width = 12, height = 7, dpi = 600)



### 1950 - 2010 ###
p <- ggplot(data = new_sf %>% drop_na(change_1950_mm_perc) %>% filter(season == 'annual'))  %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1950_mm_perc, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "Single day inc extremes \n% change 1950-2010 \n Annual", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	+ facet_wrap(. ~ return_period) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_extreme_inc_mm_perc_1950_annual.png"), p,  width = 12, height = 7, dpi = 600)



### 1970 - 2010 ###
p <- ggplot(data = new_sf %>% drop_na(change_1970_mm_perc) %>% filter(season == 'annual'))  %>%
	+ geom_sf(data = lakes, fill = "#9bbff4") %>%
	+ geom_sf(data = rivers, colour = "#4a80f5") %>%
	+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
	+ geom_sf(aes(color = change_1970_mm_perc, shape = CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "Single day inc extremes \n% change 1970-2010 \n Annual", palette = "RdBu", limits = c(-max_value, max_value), oob=squish, direction = 1) %>%
	+ scale_fill_identity() %>%
	+ facet_wrap(. ~ return_period) %>%
	+ theme_bw(12) %>%
	+ theme(legend.position = "right") %>%
	+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
	+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
	+ xlab("Longitude") %>%
	+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_extreme_inc_mm_perc_1970_annual.png"), p,  width = 12, height = 7, dpi = 600)






















