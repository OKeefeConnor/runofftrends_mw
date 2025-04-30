
# *------------------------------------------------------------------
# | FILE NAME: 02_runoff_seasonal_glm.R
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
write_figures_path <- file.path(output_path, "figures_02_runoff_seasonal_glm")
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

###########################################################################
###  increasing limb, log-flow, seasonal, GEV distribution
###########################################################################
parameterCd <- "00060"  # Discharge

### Create empty data frame
id_list <- mw_basins$STAID

model_df <- data.frame(STAID = id_list) %>%
	mutate(loc_winter = NA) %>%
	mutate(loc_spring = NA) %>%
	mutate(loc_summer = NA) %>%
	mutate(loc_autumn = NA) %>%
	mutate(loc_annual = NA) %>%
	mutate(scale_winter = NA) %>%
	mutate(scale_spring = NA) %>%
	mutate(scale_summer = NA) %>%
	mutate(scale_autumn = NA) %>%
	mutate(scale_annual = NA) %>%
	mutate(shape_winter = NA) %>%
	mutate(shape_spring = NA) %>%
	mutate(shape_summer = NA) %>%
	mutate(shape_autumn = NA) %>%
	mutate(shape_annual = NA) %>%
	mutate(p_loc_winter = NA) %>%
	mutate(p_loc_spring = NA) %>%
	mutate(p_loc_summer = NA) %>%
	mutate(p_loc_autumn = NA) %>%
	mutate(p_loc_annual = NA) %>%
	mutate(p_scale_winter = NA) %>%
	mutate(p_scale_spring = NA) %>%
	mutate(p_scale_summer = NA) %>%
	mutate(p_scale_autumn = NA) %>%
	mutate(p_scale_annual = NA) %>%
	mutate(p_shape_winter = NA) %>%
	mutate(p_shape_spring = NA) %>%
	mutate(p_shape_summer = NA) %>%
	mutate(p_shape_autumn = NA) %>%
	mutate(p_shape_annual = NA) %>%
	mutate(freq_winter = NA) %>%
	mutate(freq_spring = NA) %>%
	mutate(freq_summer = NA) %>%
	mutate(freq_autumn = NA) %>%
	mutate(freq_annual = NA) %>%
	mutate(p_freq_winter = NA) %>%
	mutate(p_freq_spring = NA) %>%
	mutate(p_freq_summer = NA) %>%
	mutate(p_freq_autumn = NA) %>%
	mutate(p_freq_annual = NA)
	
	

### Run loop for each gauge
for (k in seq(1, dim(model_df)[1])){
cat(k)
cat("\n")

	### Read in flow from USGS
	gauge_k <- model_df$STAID[k]

	### Complete flow record
	parameterCd <- "00060"  # Discharge
	startDate <- "1870-10-01"
	endDate <- "2023-09-30"

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
	  
	### Add columns for year, month, and day
	flow_k_count <- flow_k %>%
	  add_column(year = year(flow_k$date)) %>%
	  add_column(month = month(flow_k$date)) %>%
	  add_column(jdate = day(flow_k$date))


	### Filter for seasons
	flow_k_winter <- flow_k_count%>%
		filter(month == "1" | month == "2" | month == "12") %>%
		as.data.frame() %>% 
		mutate(year = case_when(month == 12 ~ year+1, 
					   TRUE  ~ year))

	flow_k_spring <- flow_k_count%>%
		filter(month == "3" | month == "4" | month == "5") %>%
		as.data.frame()

	flow_k_summer <- flow_k_count%>%
		filter(month == "6" | month == "7" | month == "8") %>%
		as.data.frame()

	flow_k_autumn <- flow_k_count%>%
		filter(month == "9" | month == "10" | month == "11") %>%
		as.data.frame()
	


	### Count number of days for each year
	flow_k_winter_count <- flow_k_winter %>%
		group_by(year) %>%
		summarize(count = n())

	flow_k_spring_count <- flow_k_spring %>%
		group_by(year) %>%
		summarize(count = n())		

	flow_k_summer_count <- flow_k_summer %>%
		group_by(year) %>%
		summarize(count = n())		

	flow_k_autumn_count <- flow_k_autumn %>%
		group_by(year) %>%
		summarize(count = n())



	### Filter to only years with more than 75 days of observations per season (no more than 2 weeks missing)
	flow_k_winter_count<- flow_k_winter_count%>% filter(count > 75)
	flow_k_spring_count <- flow_k_spring_count %>% filter(count > 75)
	flow_k_summer_count <- flow_k_summer_count %>% filter(count > 75)
	flow_k_autumn_count <- flow_k_autumn_count %>% filter(count > 75)


	### Separate the flow
	flow_sep <- process_flow_runoff(flow_k)

	### Cut to only increasing limb
	inc_df <- flow_sep$inc


	
		## Add a log10 transform
		inc_df <- inc_df %>%
			mutate(diff_log = log10(diff))


		### Filter for seasons
		inc_df_winter <- inc_df%>%
			filter(month == "1" | month == "2" | month == "12") %>%
			as.data.frame() %>% 
			mutate(year = case_when(month == 12 ~ year+1, TRUE  ~ year)) %>%
			filter(year %in% flow_k_winter_count$year)

		inc_df_spring <- inc_df%>%
			filter(month == "3" | month == "4" | month == "5") %>%
			filter(year %in% flow_k_spring_count$year) %>%
			as.data.frame()

		inc_df_summer <- inc_df%>%
			filter(month == "6" | month == "7" | month == "8") %>%
			filter(year %in% flow_k_summer_count$year) %>%
			as.data.frame()

		inc_df_autumn <- inc_df%>%
			filter(month == "9" | month == "10" | month == "11") %>%
			filter(year %in% flow_k_autumn_count$year) %>%
			as.data.frame()
		
		inc_df_annual <- inc_df %>%
			as.data.frame()



		###########################################################################
		###  Number of days with increasing flow (*Not used in "Trends in Seasonal Streamflow and Daily Runoff... "*)
		###########################################################################
		first_year <- min(inc_df$year, na.rm=TRUE)
		last_year <- max(inc_df$year, na.rm=TRUE)
		
		n_inc_df_winter = NULL
		n_inc_df_spring = NULL
		n_inc_df_summer = NULL
		n_inc_df_autumn = NULL
		n_inc_df_annual = NULL

		for (x in seq(first_year, last_year)) { #Run for loop for each year
			year = x
			
			### count number of days with increasing flow in each season
			n_inc_winter = length(which(inc_df_winter$year == x))
			n_inc_spring = length(which(inc_df_spring$year == x))
			n_inc_summer = length(which(inc_df_summer$year == x))
			n_inc_autumn = length(which(inc_df_autumn$year == x))
			n_inc_annual = length(which(inc_df_annual$year == x))
			
			### create dataframe with column for year and column for number of days 	(filter out years without data)
			n_inc_df_winter = rbind(n_inc_df_winter, data.frame(year,n_inc_winter)) %>% filter(n_inc_winter > 0)
			n_inc_df_spring = rbind(n_inc_df_spring, data.frame(year,n_inc_spring)) %>% filter(n_inc_spring > 0)
			n_inc_df_summer = rbind(n_inc_df_summer, data.frame(year,n_inc_summer)) %>% filter(n_inc_summer > 0)
			n_inc_df_autumn = rbind(n_inc_df_autumn, data.frame(year,n_inc_autumn)) %>% filter(n_inc_autumn > 0)
			n_inc_df_annual = rbind(n_inc_df_annual, data.frame(year,n_inc_annual)) %>% filter(n_inc_annual > 0)

		}
		

		
		### Linear regression for number of days with increasing flow for each season
		lm_winter_n_inc <- lm(n_inc_winter ~ year, n_inc_df_winter)
		freq_winter <- summary(lm_winter_n_inc)$coefficients["year", "Estimate"]
		p_freq_winter <- summary(lm_winter_n_inc)$coefficients[2,4]  
		
		lm_spring_n_inc <- lm(n_inc_spring ~ year, n_inc_df_spring)
		freq_spring <- summary(lm_spring_n_inc)$coefficients["year", "Estimate"]
		p_freq_spring <- summary(lm_spring_n_inc)$coefficients[2,4]  
		
		lm_summer_n_inc <- lm(n_inc_summer ~ year, n_inc_df_summer)
		freq_summer <- summary(lm_summer_n_inc)$coefficients["year", "Estimate"]
		p_freq_summer <- summary(lm_summer_n_inc)$coefficients[2,4]  
		
		lm_autumn_n_inc <- lm(n_inc_autumn ~ year, n_inc_df_autumn)
		freq_autumn <- summary(lm_autumn_n_inc)$coefficients["year", "Estimate"]
		p_freq_autumn <- summary(lm_autumn_n_inc)$coefficients[2,4]  
		
		lm_annual_n_inc <- lm(n_inc_annual ~ year, n_inc_df_annual)
		freq_annual <- summary(lm_annual_n_inc)$coefficients["year", "Estimate"]
		p_freq_annual <- summary(lm_annual_n_inc)$coefficients[2,4]  


		### Save into columns
		model_df$freq_winter[k] <- freq_winter
		model_df$freq_spring[k] <- freq_spring
		model_df$freq_summer[k] <- freq_summer
		model_df$freq_autumn[k] <- freq_autumn
		model_df$freq_annual[k] <- freq_annual
		
		model_df$p_freq_winter[k] <- p_freq_winter
		model_df$p_freq_spring[k] <- p_freq_spring
		model_df$p_freq_summer[k] <- p_freq_summer
		model_df$p_freq_autumn[k] <- p_freq_autumn
		model_df$p_freq_annual[k] <- p_freq_annual


		###########################################################################
		###  Fit a GEV for each season
		###########################################################################

		### Fit gevlss model for each season
		gev_fit_winter <- gam(list(diff_log 
		 ~ year,
		 ~ year,
		 ~ year),
		family=gevlss(link = list("identity", "identity", "identity")), 
		data = inc_df_winter,
		select = TRUE,
		optimizer="efs")

		gev_fit_spring <- gam(list(diff_log
		 ~ year,
		 ~ year,
		 ~ year),
		family=gevlss(link = list("identity", "identity", "identity")), 
		data = inc_df_spring,
		select = TRUE,
		optimizer="efs")

		gev_fit_summer <- gam(list(diff_log
		 ~ year,
		 ~ year,
		 ~ year),			
		family=gevlss(link = list("identity", "identity", "identity")), 
		data = inc_df_summer,
		select = TRUE,
		optimizer="efs")

		gev_fit_autumn <- gam(list(diff_log
		 ~ year,
		 ~ year,
		 ~ year),			
		family=gevlss(link = list("identity", "identity", "identity")), 
		data = inc_df_autumn,
		select = TRUE,
		optimizer="efs")

		gev_fit_annual <- gam(list(diff_log
		 ~ year,
		 ~ year,
		 ~ year),			
		family=gevlss(link = list("identity", "identity", "identity")), 
		data = inc_df_annual,
		select = TRUE,
		optimizer="efs")
  


		###########################################################################
		###  Generate a full estimate
		###########################################################################
		
		### Create dataframe to find parameters for every 10 years
		###year_df <- data.frame(year = seq(min(flow_k_autumn_avg$year), max(flow_k_autumn_avg$year), by = 10))
		year_df <- data.frame(year = seq(1900, 2020, by = 10))

		### Randomly draw from the estimates
		gev_draws_winter <- draw_from_gev_linear(mod = gev_fit_winter, newdata = year_df, n = 1000)

		gev_draws_spring <- draw_from_gev_linear(mod = gev_fit_spring, newdata = year_df, n = 1000)
		
		gev_draws_summer <- draw_from_gev_linear(mod = gev_fit_summer, newdata = year_df, n = 1000)

		gev_draws_autumn <- draw_from_gev_linear(mod = gev_fit_autumn, newdata = year_df, n = 1000)
		
		gev_draws_annual <- draw_from_gev_linear(mod = gev_fit_annual, newdata = year_df, n = 1000)


		### Extract an estimate of the scale parameter
		scale_est_winter <- summarize_gev_draw(gev_draws_winter, var = "scale", conf_int = 0.95)
		predict_scale_winter_temp <- scale_est_winter$summary
		predict_scale_winter_temp <- predict_scale_winter_temp %>%
			mutate(scale_winter_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(scale_winter = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper) %>%
			mutate(STAID = gauge_k)

		scale_est_spring <- summarize_gev_draw(gev_draws_spring, var = "scale", conf_int = 0.95)
		predict_scale_spring_temp <- scale_est_spring$summary
		predict_scale_spring_temp <- predict_scale_spring_temp %>%
			mutate(scale_spring_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(scale_spring = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper)
			
		scale_est_summer <- summarize_gev_draw(gev_draws_summer, var = "scale", conf_int = 0.95)
		predict_scale_summer_temp <- scale_est_summer$summary
		predict_scale_summer_temp <- predict_scale_summer_temp %>%
			mutate(scale_summer_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(scale_summer = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper)

		scale_est_autumn <- summarize_gev_draw(gev_draws_autumn, var = "scale", conf_int = 0.95)
		predict_scale_autumn_temp <- scale_est_autumn$summary
		predict_scale_autumn_temp <- predict_scale_autumn_temp %>%
			mutate(scale_autumn_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(scale_autumn = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper)
	
		scale_est_annual <- summarize_gev_draw(gev_draws_annual, var = "scale", conf_int = 0.95)
		predict_scale_annual_temp <- scale_est_annual$summary
		predict_scale_annual_temp <- predict_scale_annual_temp %>%
			mutate(scale_annual_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(scale_annual = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper)



		### Extract an estimate of the shape parameter
		shape_est_winter <- summarize_gev_draw(gev_draws_winter, var = "shape", conf_int = 0.95)
		predict_shape_winter_temp <- shape_est_winter$summary
		predict_shape_winter_temp <- predict_shape_winter_temp %>%
			mutate(shape_winter_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(shape_winter = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper)

		shape_est_spring <- summarize_gev_draw(gev_draws_spring, var = "shape", conf_int = 0.95)
		predict_shape_spring_temp <- shape_est_spring$summary
		predict_shape_spring_temp <- predict_shape_spring_temp %>%
			mutate(shape_spring_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(shape_spring = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper)

		shape_est_summer <- summarize_gev_draw(gev_draws_summer, var = "shape", conf_int = 0.95)
		predict_shape_summer_temp <- shape_est_summer$summary
		predict_shape_summer_temp <- predict_shape_summer_temp %>%
			mutate(shape_summer_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(shape_summer = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper)

		shape_est_autumn <- summarize_gev_draw(gev_draws_autumn, var = "shape", conf_int = 0.95)
		predict_shape_autumn_temp <- shape_est_autumn$summary
		predict_shape_autumn_temp <- predict_shape_autumn_temp %>%
			mutate(shape_autumn_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(shape_autumn = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper)

		shape_est_annual <- summarize_gev_draw(gev_draws_annual, var = "shape", conf_int = 0.95)
		predict_shape_annual_temp <- shape_est_annual$summary
		predict_shape_annual_temp <- predict_shape_annual_temp %>%
			mutate(shape_annual_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(shape_annual = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper)

 
 
 		### Extract an estimate of the location parameter
		loc_est_winter <- summarize_gev_draw(gev_draws_winter, var = "loc", conf_int = 0.95)
		predict_loc_winter_temp <- loc_est_winter$summary
		predict_loc_winter_temp <- predict_loc_winter_temp %>%
			mutate(loc_winter_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(loc_winter = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper)

		loc_est_spring <- summarize_gev_draw(gev_draws_spring, var = "loc", conf_int = 0.95)
		predict_loc_spring_temp <- loc_est_spring$summary
		predict_loc_spring_temp <- predict_loc_spring_temp %>%
			mutate(loc_spring_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(loc_spring = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper)

		loc_est_summer <- summarize_gev_draw(gev_draws_summer, var = "loc", conf_int = 0.95)
		predict_loc_summer_temp <- loc_est_summer$summary
		predict_loc_summer_temp <- predict_loc_summer_temp %>%
			mutate(loc_summer_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(loc_summer = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper)

		loc_est_autumn <- summarize_gev_draw(gev_draws_autumn, var = "loc", conf_int = 0.95)
		predict_loc_autumn_temp <- loc_est_autumn$summary
		predict_loc_autumn_temp <- predict_loc_autumn_temp %>%
			mutate(loc_autumn_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(loc_autumn = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper)

		loc_est_annual <- summarize_gev_draw(gev_draws_annual, var = "loc", conf_int = 0.95)
		predict_loc_annual_temp <- loc_est_annual$summary
		predict_loc_annual_temp <- predict_loc_annual_temp %>%
			mutate(loc_annual_std = (est_median-est_median[8])/est_median[8]) %>%
			rename(loc_annual = est_median) %>%
			select(-est_lower,-est_upper,-slope_median,-slope_lower,-slope_upper)



		### Extract distribution estimates for quantiles (25 = low flow, 75 = high flow, 50 = median flow)
		dist_low_winter <- summarize_gev_draw(gev_draws_winter, var = "quantile", level = 0.25, conf_int = 0.95)$summary
		dist_high_winter <- summarize_gev_draw(gev_draws_winter, var = "quantile", level = 0.75, conf_int = 0.95)$summary
		dist_med_winter <- summarize_gev_draw(gev_draws_winter, var = "quantile", level = 0.50, conf_int = 0.95)$summary

		dist_low_spring <- summarize_gev_draw(gev_draws_spring, var = "quantile", level = 0.25, conf_int = 0.95)$summary
		dist_high_spring <- summarize_gev_draw(gev_draws_spring, var = "quantile", level = 0.75, conf_int = 0.95)$summary
		dist_med_spring <- summarize_gev_draw(gev_draws_spring, var = "quantile", level = 0.50, conf_int = 0.95)$summary

		dist_low_summer <- summarize_gev_draw(gev_draws_summer, var = "quantile", level = 0.25, conf_int = 0.95)$summary
		dist_high_summer <- summarize_gev_draw(gev_draws_summer, var = "quantile", level = 0.75, conf_int = 0.95)$summary
		dist_med_summer <- summarize_gev_draw(gev_draws_summer, var = "quantile", level = 0.50, conf_int = 0.95)$summary

		dist_low_autumn <- summarize_gev_draw(gev_draws_autumn, var = "quantile", level = 0.25, conf_int = 0.95)$summary
		dist_high_autumn <- summarize_gev_draw(gev_draws_autumn, var = "quantile", level = 0.75, conf_int = 0.95)$summary
		dist_med_autumn <- summarize_gev_draw(gev_draws_autumn, var = "quantile", level = 0.50, conf_int = 0.95)$summary

		dist_low_annual <- summarize_gev_draw(gev_draws_annual, var = "quantile", level = 0.25, conf_int = 0.95)$summary
		dist_high_annual <- summarize_gev_draw(gev_draws_annual, var = "quantile", level = 0.75, conf_int = 0.95)$summary
		dist_med_annual <- summarize_gev_draw(gev_draws_annual, var = "quantile", level = 0.50, conf_int = 0.95)$summary

		quantile_winter_df_temp <- data.frame(year = dist_low_winter$year, low_25_winter = dist_low_winter$est_median, high_75_winter = dist_high_winter$est_median, med_50_winter = dist_med_winter$est_median,
			low_slope_median = dist_low_winter$slope_median, high_slope_median = dist_high_winter$slope_median, med_slope_median = dist_med_winter$slope_median,
			low_slope_upper = dist_low_winter$slope_upper, high_slope_upper = dist_high_winter$slope_upper, med_slope_upper = dist_med_winter$slope_upper,
			low_slope_lower = dist_low_winter$slope_lower, high_slope_lower = dist_high_winter$slope_lower, med_slope_lower = dist_med_winter$slope_lower)
		quantile_winter_df_temp <- quantile_winter_df_temp %>%
			mutate(low_25_winter_std = (low_25_winter-low_25_winter[8])/low_25_winter[8]) %>%
			mutate(high_75_winter_std = (high_75_winter-high_75_winter[8])/high_75_winter[8]) %>%
			mutate(med_50_winter_std = (med_50_winter-med_50_winter[8])/med_50_winter[8]) %>%
			mutate(low_25_winter_slope_median = low_slope_median[8]) %>%
			mutate(low_25_winter_slope_upper = low_slope_upper[8]) %>%
			mutate(low_25_winter_slope_lower = low_slope_lower[8]) %>%
			mutate(high_75_winter_slope_median = high_slope_median[8]) %>%
			mutate(high_75_winter_slope_upper = high_slope_upper[8]) %>%
			mutate(high_75_winter_slope_lower = high_slope_lower[8]) %>%
			mutate(med_50_winter_slope_median = med_slope_median[8]) %>%
			mutate(med_50_winter_slope_upper = med_slope_upper[8]) %>%
			mutate(med_50_winter_slope_lower = med_slope_lower[8]) %>%
			select(-low_25_winter,-high_75_winter,-med_50_winter,-low_slope_median,-high_slope_median,-med_slope_median,-low_slope_upper,-high_slope_upper,-med_slope_upper,-low_slope_lower,-high_slope_lower,-med_slope_lower)

		quantile_spring_df_temp <- data.frame(year = dist_low_spring$year, low_25_spring = dist_low_spring$est_median, high_75_spring = dist_high_spring$est_median, med_50_spring = dist_med_spring$est_median,
			low_slope_median = dist_low_spring$slope_median, high_slope_median = dist_high_spring$slope_median, med_slope_median = dist_med_spring$slope_median,
			low_slope_upper = dist_low_spring$slope_upper, high_slope_upper = dist_high_spring$slope_upper, med_slope_upper = dist_med_spring$slope_upper,
			low_slope_lower = dist_low_spring$slope_lower, high_slope_lower = dist_high_spring$slope_lower, med_slope_lower = dist_med_spring$slope_lower)
		quantile_spring_df_temp <- quantile_spring_df_temp %>%
			mutate(low_25_spring_std = (low_25_spring-low_25_spring[8])/low_25_spring[8]) %>%
			mutate(high_75_spring_std = (high_75_spring-high_75_spring[8])/high_75_spring[8]) %>%
			mutate(med_50_spring_std = (med_50_spring-med_50_spring[8])/med_50_spring[8]) %>%
			mutate(low_25_spring_slope_median = low_slope_median[8]) %>%
			mutate(low_25_spring_slope_upper = low_slope_upper[8]) %>%
			mutate(low_25_spring_slope_lower = low_slope_lower[8]) %>%
			mutate(high_75_spring_slope_median = high_slope_median[8]) %>%
			mutate(high_75_spring_slope_upper = high_slope_upper[8]) %>%
			mutate(high_75_spring_slope_lower = high_slope_lower[8]) %>%
			mutate(med_50_spring_slope_median = med_slope_median[8]) %>%
			mutate(med_50_spring_slope_upper = med_slope_upper[8]) %>%
			mutate(med_50_spring_slope_lower = med_slope_lower[8]) %>%
			select(-low_25_spring,-high_75_spring,-med_50_spring,-low_slope_median,-high_slope_median,-med_slope_median,-low_slope_upper,-high_slope_upper,-med_slope_upper,-low_slope_lower,-high_slope_lower,-med_slope_lower)

		quantile_summer_df_temp <- data.frame(year = dist_low_summer$year, low_25_summer = dist_low_summer$est_median, high_75_summer = dist_high_summer$est_median, med_50_summer = dist_med_summer$est_median,
			low_slope_median = dist_low_summer$slope_median, high_slope_median = dist_high_summer$slope_median, med_slope_median = dist_med_summer$slope_median,
			low_slope_upper = dist_low_summer$slope_upper, high_slope_upper = dist_high_summer$slope_upper, med_slope_upper = dist_med_summer$slope_upper,
			low_slope_lower = dist_low_summer$slope_lower, high_slope_lower = dist_high_summer$slope_lower, med_slope_lower = dist_med_summer$slope_lower)
		quantile_summer_df_temp <- quantile_summer_df_temp %>%
			mutate(low_25_summer_std = (low_25_summer-low_25_summer[8])/low_25_summer[8]) %>%
			mutate(high_75_summer_std = (high_75_summer-high_75_summer[8])/high_75_summer[8]) %>%
			mutate(med_50_summer_std = (med_50_summer-med_50_summer[8])/med_50_summer[8]) %>%
			mutate(low_25_summer_slope_median = low_slope_median[8]) %>%
			mutate(low_25_summer_slope_upper = low_slope_upper[8]) %>%
			mutate(low_25_summer_slope_lower = low_slope_lower[8]) %>%
			mutate(high_75_summer_slope_median = high_slope_median[8]) %>%
			mutate(high_75_summer_slope_upper = high_slope_upper[8]) %>%
			mutate(high_75_summer_slope_lower = high_slope_lower[8]) %>%
			mutate(med_50_summer_slope_median = med_slope_median[8]) %>%
			mutate(med_50_summer_slope_upper = med_slope_upper[8]) %>%
			mutate(med_50_summer_slope_lower = med_slope_lower[8]) %>%
			select(-low_25_summer,-high_75_summer,-med_50_summer,-low_slope_median,-high_slope_median,-med_slope_median,-low_slope_upper,-high_slope_upper,-med_slope_upper,-low_slope_lower,-high_slope_lower,-med_slope_lower)

		quantile_autumn_df_temp <- data.frame(year = dist_low_autumn$year, low_25_autumn = dist_low_autumn$est_median, high_75_autumn = dist_high_autumn$est_median, med_50_autumn = dist_med_autumn$est_median,
			low_slope_median = dist_low_autumn$slope_median, high_slope_median = dist_high_autumn$slope_median, med_slope_median = dist_med_autumn$slope_median,
			low_slope_upper = dist_low_autumn$slope_upper, high_slope_upper = dist_high_autumn$slope_upper, med_slope_upper = dist_med_autumn$slope_upper,
			low_slope_lower = dist_low_autumn$slope_lower, high_slope_lower = dist_high_autumn$slope_lower, med_slope_lower = dist_med_autumn$slope_lower)
		quantile_autumn_df_temp <- quantile_autumn_df_temp %>%
			mutate(low_25_autumn_std = (low_25_autumn-low_25_autumn[8])/low_25_autumn[8]) %>%
			mutate(high_75_autumn_std = (high_75_autumn-high_75_autumn[8])/high_75_autumn[8]) %>%
			mutate(med_50_autumn_std = (med_50_autumn-med_50_autumn[8])/med_50_autumn[8]) %>%
			mutate(low_25_autumn_slope_median = low_slope_median[8]) %>%
			mutate(low_25_autumn_slope_upper = low_slope_upper[8]) %>%
			mutate(low_25_autumn_slope_lower = low_slope_lower[8]) %>%
			mutate(high_75_autumn_slope_median = high_slope_median[8]) %>%
			mutate(high_75_autumn_slope_upper = high_slope_upper[8]) %>%
			mutate(high_75_autumn_slope_lower = high_slope_lower[8]) %>%
			mutate(med_50_autumn_slope_median = med_slope_median[8]) %>%
			mutate(med_50_autumn_slope_upper = med_slope_upper[8]) %>%
			mutate(med_50_autumn_slope_lower = med_slope_lower[8]) %>%
			select(-low_25_autumn,-high_75_autumn,-med_50_autumn,-low_slope_median,-high_slope_median,-med_slope_median,-low_slope_upper,-high_slope_upper,-med_slope_upper,-low_slope_lower,-high_slope_lower,-med_slope_lower)

		quantile_annual_df_temp <- data.frame(year = dist_low_annual$year, low_25_annual = dist_low_annual$est_median, high_75_annual = dist_high_annual$est_median, med_50_annual = dist_med_annual$est_median,
			low_slope_median = dist_low_annual$slope_median, high_slope_median = dist_high_annual$slope_median, med_slope_median = dist_med_annual$slope_median,
			low_slope_upper = dist_low_annual$slope_upper, high_slope_upper = dist_high_annual$slope_upper, med_slope_upper = dist_med_annual$slope_upper,
			low_slope_lower = dist_low_annual$slope_lower, high_slope_lower = dist_high_annual$slope_lower, med_slope_lower = dist_med_annual$slope_lower)
		quantile_annual_df_temp <- quantile_annual_df_temp %>%
			mutate(low_25_annual_std = (low_25_annual-low_25_annual[8])/low_25_annual[8]) %>%
			mutate(high_75_annual_std = (high_75_annual-high_75_annual[8])/high_75_annual[8]) %>%
			mutate(med_50_annual_std = (med_50_annual-med_50_annual[8])/med_50_annual[8]) %>%
			mutate(low_25_annual_slope_median = low_slope_median[8]) %>%
			mutate(low_25_annual_slope_upper = low_slope_upper[8]) %>%
			mutate(low_25_annual_slope_lower = low_slope_lower[8]) %>%
			mutate(high_75_annual_slope_median = high_slope_median[8]) %>%
			mutate(high_75_annual_slope_upper = high_slope_upper[8]) %>%
			mutate(high_75_annual_slope_lower = high_slope_lower[8]) %>%
			mutate(med_50_annual_slope_median = med_slope_median[8]) %>%
			mutate(med_50_annual_slope_upper = med_slope_upper[8]) %>%
			mutate(med_50_annual_slope_lower = med_slope_lower[8]) %>%
			select(-low_25_annual,-high_75_annual,-med_50_annual,-low_slope_median,-high_slope_median,-med_slope_median,-low_slope_upper,-high_slope_upper,-med_slope_upper,-low_slope_lower,-high_slope_lower,-med_slope_lower)

		
		
		### Determine first and last year
		first_year <- min(inc_df$year, na.rm=TRUE)
		last_year <- max(inc_df$year, na.rm=TRUE)
		#record_length <- last_year - first_year + 1
		
		
		### Create temporary dataframe
		plot_df_temp <- predict_scale_winter_temp		
		## Filter dataframe based on years of avaialable data
		plot_df_temp <- plot_df_temp %>% filter(year >= first_year & year <= last_year)		
		## Join the tempoerary dataframes together
		plot_df_temp <- plot_df_temp %>% left_join(predict_scale_spring_temp, by = 'year')
		plot_df_temp <- plot_df_temp %>% left_join(predict_scale_summer_temp, by = 'year')
		plot_df_temp <- plot_df_temp %>% left_join(predict_scale_autumn_temp, by = 'year')
		plot_df_temp <- plot_df_temp %>% left_join(predict_scale_annual_temp, by = 'year')
		
		plot_df_temp <- plot_df_temp %>% left_join(predict_shape_winter_temp, by = 'year')
		plot_df_temp <- plot_df_temp %>% left_join(predict_shape_spring_temp, by = 'year')
		plot_df_temp <- plot_df_temp %>% left_join(predict_shape_summer_temp, by = 'year')
		plot_df_temp <- plot_df_temp %>% left_join(predict_shape_autumn_temp, by = 'year')
		plot_df_temp <- plot_df_temp %>% left_join(predict_shape_annual_temp, by = 'year')
		
		plot_df_temp <- plot_df_temp %>% left_join(predict_loc_winter_temp, by = 'year')
		plot_df_temp <- plot_df_temp %>% left_join(predict_loc_spring_temp, by = 'year')
		plot_df_temp <- plot_df_temp %>% left_join(predict_loc_summer_temp, by = 'year')
		plot_df_temp <- plot_df_temp %>% left_join(predict_loc_autumn_temp, by = 'year')
		plot_df_temp <- plot_df_temp %>% left_join(predict_loc_annual_temp, by = 'year')
		
		plot_df_temp <- plot_df_temp %>% left_join(quantile_winter_df_temp, by = 'year')
		plot_df_temp <- plot_df_temp %>% left_join(quantile_spring_df_temp, by = 'year')
		plot_df_temp <- plot_df_temp %>% left_join(quantile_summer_df_temp, by = 'year')
		plot_df_temp <- plot_df_temp %>% left_join(quantile_autumn_df_temp, by = 'year')	
		plot_df_temp <- plot_df_temp %>% left_join(quantile_annual_df_temp, by = 'year')



		### Extract parameters and p-values
		loc_winter <- summary(gev_fit_winter)$p.table[2, 1]
		scale_winter <- summary(gev_fit_winter)$p.table[4, 1]
		shape_winter <- summary(gev_fit_winter)$p.table[6, 1]
		p_loc_winter <- summary(gev_fit_winter)$p.table[2, 4]
		p_scale_winter <- summary(gev_fit_winter)$p.table[4, 4]
		p_shape_winter <- summary(gev_fit_winter)$p.table[6, 4]
		
		loc_spring <- summary(gev_fit_spring)$p.table[2, 1]
		scale_spring <- summary(gev_fit_spring)$p.table[4, 1]
		shape_spring <- summary(gev_fit_spring)$p.table[6, 1]
		p_loc_spring <- summary(gev_fit_spring)$p.table[2, 4]
		p_scale_spring <- summary(gev_fit_spring)$p.table[4, 4]
		p_shape_spring <- summary(gev_fit_spring)$p.table[6, 4]

		loc_summer <- summary(gev_fit_summer)$p.table[2, 1]
		scale_summer <- summary(gev_fit_summer)$p.table[4, 1]
		shape_summer <- summary(gev_fit_summer)$p.table[6, 1]
		p_loc_summer <- summary(gev_fit_summer)$p.table[2, 4]
		p_scale_summer <- summary(gev_fit_summer)$p.table[4, 4]
		p_shape_summer <- summary(gev_fit_summer)$p.table[6, 4]
		
		loc_autumn <- summary(gev_fit_autumn)$p.table[2, 1]
		scale_autumn <- summary(gev_fit_autumn)$p.table[4, 1]
		shape_autumn <- summary(gev_fit_autumn)$p.table[6, 1]
		p_loc_autumn <- summary(gev_fit_autumn)$p.table[2, 4]
		p_scale_autumn <- summary(gev_fit_autumn)$p.table[4, 4]
		p_shape_autumn <- summary(gev_fit_autumn)$p.table[6, 4]
		
		loc_annual <- summary(gev_fit_annual)$p.table[2, 1]
		scale_annual <- summary(gev_fit_annual)$p.table[4, 1]
		shape_annual <- summary(gev_fit_annual)$p.table[6, 1]
		p_loc_annual <- summary(gev_fit_annual)$p.table[2, 4]
		p_scale_annual <- summary(gev_fit_annual)$p.table[4, 4]
		p_shape_annual <- summary(gev_fit_annual)$p.table[6, 4]
		
		### Save into columns
		model_df$loc_winter[k] <- loc_winter
		model_df$scale_winter[k] <- scale_winter
		model_df$shape_winter[k] <- shape_winter
		model_df$p_loc_winter[k] <- p_loc_winter
		model_df$p_scale_winter[k] <- p_scale_winter
		model_df$p_shape_winter[k] <- p_shape_winter
		
		model_df$loc_spring[k] <- loc_spring
		model_df$scale_spring[k] <- scale_spring
		model_df$shape_spring[k] <- shape_spring
		model_df$p_loc_spring[k] <- p_loc_spring
		model_df$p_scale_spring[k] <- p_scale_spring
		model_df$p_shape_spring[k] <- p_shape_spring

		model_df$loc_summer[k] <- loc_summer
		model_df$scale_summer[k] <- scale_summer
		model_df$shape_summer[k] <- shape_summer
		model_df$p_loc_summer[k] <- p_loc_summer
		model_df$p_scale_summer[k] <- p_scale_summer
		model_df$p_shape_summer[k] <- p_shape_summer

		model_df$loc_autumn[k] <- loc_autumn
		model_df$scale_autumn[k] <- scale_autumn
		model_df$shape_autumn[k] <- shape_autumn
		model_df$p_loc_autumn[k] <- p_loc_autumn
		model_df$p_scale_autumn[k] <- p_scale_autumn
		model_df$p_shape_autumn[k] <- p_shape_autumn
		
		model_df$loc_annual[k] <- loc_annual
		model_df$scale_annual[k] <- scale_annual
		model_df$shape_annual[k] <- shape_annual
		model_df$p_loc_annual[k] <- p_loc_annual
		model_df$p_scale_annual[k] <- p_scale_annual
		model_df$p_shape_annual[k] <- p_shape_annual

		### Save columns into plot_df
		plot_df_temp$p_loc_winter <- p_loc_winter
		plot_df_temp$p_scale_winter <- p_scale_winter
		plot_df_temp$p_shape_winter <- p_shape_winter
		
		plot_df_temp$p_loc_spring <- p_loc_spring
		plot_df_temp$p_scale_spring <- p_scale_spring
		plot_df_temp$p_shape_spring <- p_shape_spring

		plot_df_temp$p_loc_summer <- p_loc_summer
		plot_df_temp$p_scale_summer <- p_scale_summer
		plot_df_temp$p_shape_summer <- p_shape_summer

		plot_df_temp$p_loc_autumn <- p_loc_autumn
		plot_df_temp$p_scale_autumn <- p_scale_autumn
		plot_df_temp$p_shape_autumn <- p_shape_autumn

		plot_df_temp$p_loc_annual <- p_loc_annual
		plot_df_temp$p_scale_annual <- p_scale_annual
		plot_df_temp$p_shape_annual <- p_shape_annual



if (k == 1){	#first pass through loop
	plot_df <- plot_df_temp
	} else {
		plot_df <- plot_df %>%
			bind_rows(plot_df_temp)
	}

}



### Add columns that do if then statements based on change and p-value
model_df <- model_df %>%
	mutate(change_winter_loc = NA) %>%
	mutate(change_spring_loc = NA) %>%
	mutate(change_summer_loc = NA) %>%
	mutate(change_autumn_loc = NA) %>%
	mutate(change_annual_loc = NA) %>%
	mutate(change_winter_loc = case_when(
  	loc_winter > 0 & p_loc_winter < 0.05 ~ "Increase",
  	loc_winter < 0 & p_loc_winter < 0.05 ~ "Decrease",
  	p_loc_winter >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_spring_loc = case_when(
  	loc_spring > 0 & p_loc_spring < 0.05 ~ "Increase",
  	loc_spring < 0 & p_loc_spring < 0.05 ~ "Decrease",
  	p_loc_spring >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_summer_loc = case_when(
  	loc_summer > 0 & p_loc_summer < 0.05 ~ "Increase",
  	loc_summer < 0 & p_loc_summer < 0.05 ~ "Decrease",
  	p_loc_summer >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_autumn_loc = case_when(
  	loc_autumn > 0 & p_loc_autumn < 0.05 ~ "Increase",
  	loc_autumn < 0 & p_loc_autumn < 0.05 ~ "Decrease",
  	p_loc_autumn >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_annual_loc = case_when(
  	loc_annual > 0 & p_loc_annual < 0.05 ~ "Increase",
  	loc_annual < 0 & p_loc_annual < 0.05 ~ "Decrease",
  	p_loc_annual >= 0.05 ~ "NotSig",
  	.default = NA_character_
	))

model_df <- model_df %>%
	mutate(change_winter_scale = NA) %>%
	mutate(change_spring_scale = NA) %>%
	mutate(change_summer_scale = NA) %>%
	mutate(change_autumn_scale = NA) %>%
	mutate(change_annual_scale = NA) %>%
	mutate(change_winter_scale = case_when(
  	scale_winter > 0 & p_scale_winter < 0.05 ~ "Increase",
  	scale_winter < 0 & p_scale_winter < 0.05 ~ "Decrease",
  	p_scale_winter >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_spring_scale = case_when(
  	scale_spring > 0 & p_scale_spring < 0.05 ~ "Increase",
  	scale_spring < 0 & p_scale_spring < 0.05 ~ "Decrease",
  	p_scale_spring >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_summer_scale = case_when(
  	scale_summer > 0 & p_scale_summer < 0.05 ~ "Increase",
  	scale_summer < 0 & p_scale_summer < 0.05 ~ "Decrease",
  	p_scale_summer >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_autumn_scale = case_when(
  	scale_autumn > 0 & p_scale_autumn < 0.05 ~ "Increase",
  	scale_autumn < 0 & p_scale_autumn < 0.05 ~ "Decrease",
  	p_scale_autumn >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_annual_scale = case_when(
  	scale_annual > 0 & p_scale_annual < 0.05 ~ "Increase",
  	scale_annual < 0 & p_scale_annual < 0.05 ~ "Decrease",
  	p_scale_annual >= 0.05 ~ "NotSig",
  	.default = NA_character_
	))

model_df <- model_df %>%
	mutate(change_winter_shape = NA) %>%
	mutate(change_spring_shape = NA) %>%
	mutate(change_summer_shape = NA) %>%
	mutate(change_autumn_shape = NA) %>%
	mutate(change_annual_shape = NA) %>%
	mutate(change_winter_shape = case_when(
  	shape_winter > 0 & p_shape_winter < 0.05 ~ "Increase",
  	shape_winter < 0 & p_shape_winter < 0.05 ~ "Decrease",
  	p_shape_winter >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_spring_shape = case_when(
  	shape_spring > 0 & p_shape_spring < 0.05 ~ "Increase",
  	shape_spring < 0 & p_shape_spring < 0.05 ~ "Decrease",
  	p_shape_spring >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_summer_shape = case_when(
  	shape_summer > 0 & p_shape_summer < 0.05 ~ "Increase",
  	shape_summer < 0 & p_shape_summer < 0.05 ~ "Decrease",
  	p_shape_summer >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_autumn_shape = case_when(
  	shape_autumn > 0 & p_shape_autumn < 0.05 ~ "Increase",
  	shape_autumn < 0 & p_shape_autumn < 0.05 ~ "Decrease",
  	p_shape_autumn >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_annual_shape = case_when(
  	shape_annual > 0 & p_shape_annual < 0.05 ~ "Increase",
  	shape_annual < 0 & p_shape_annual < 0.05 ~ "Decrease",
  	p_shape_annual >= 0.05 ~ "NotSig",
  	.default = NA_character_
	))

model_df <- model_df %>%
	mutate(change_winter_freq = NA) %>%
	mutate(change_spring_freq = NA) %>%
	mutate(change_summer_freq = NA) %>%
	mutate(change_autumn_freq = NA) %>%
	mutate(change_annual_freq = NA) %>%
	mutate(change_winter_freq = case_when(
  	freq_winter > 0 & p_freq_winter < 0.05 ~ "Increase",
  	freq_winter < 0 & p_freq_winter < 0.05 ~ "Decrease",
  	p_freq_winter >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_spring_freq = case_when(
  	freq_spring > 0 & p_freq_spring < 0.05 ~ "Increase",
  	freq_spring < 0 & p_freq_spring < 0.05 ~ "Decrease",
  	p_freq_spring >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_summer_freq = case_when(
  	freq_summer > 0 & p_freq_summer < 0.05 ~ "Increase",
  	freq_summer < 0 & p_freq_summer < 0.05 ~ "Decrease",
  	p_freq_summer >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_autumn_freq = case_when(
  	freq_autumn > 0 & p_freq_autumn < 0.05 ~ "Increase",
  	freq_autumn < 0 & p_freq_autumn < 0.05 ~ "Decrease",
  	p_freq_autumn >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_annual_freq = case_when(
  	freq_annual > 0 & p_freq_annual < 0.05 ~ "Increase",
  	freq_annual < 0 & p_freq_annual < 0.05 ~ "Decrease",
  	p_freq_annual >= 0.05 ~ "NotSig",
  	.default = NA_character_
	))



plot_df <- plot_df %>%
	mutate(change_winter_75 = NA) %>%
	mutate(change_winter_25 = NA) %>%
	mutate(change_winter_50 = NA) %>%
	mutate(change_winter_IQR = NA) %>%
	mutate(change_spring_75 = NA) %>%
	mutate(change_spring_25 = NA) %>%
	mutate(change_spring_50 = NA) %>%
	mutate(change_spring_IQR = NA) %>%
	mutate(change_summer_75 = NA) %>%
	mutate(change_summer_25 = NA) %>%
	mutate(change_summer_50 = NA) %>%
	mutate(change_summer_IQR = NA) %>%
	mutate(change_autumn_75 = NA) %>%
	mutate(change_autumn_25 = NA) %>%
	mutate(change_autumn_50 = NA) %>%
	mutate(change_autumn_IQR = NA) %>%
	mutate(change_annual_75 = NA) %>%
	mutate(change_annual_25 = NA) %>%
	mutate(change_annual_50 = NA) %>%
	mutate(change_annual_IQR = NA) %>%
	
	mutate(change_winter_75 = case_when(
  	high_75_winter_slope_median > 0 & high_75_winter_slope_lower > 0 ~ "Increase",
  	high_75_winter_slope_median < 0 & high_75_winter_slope_upper < 0 ~ "Decrease",
  	high_75_winter_slope_median > 0 & high_75_winter_slope_lower < 0 ~ "NotSig",
	high_75_winter_slope_median < 0 & high_75_winter_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_winter_25 = case_when(
  	low_25_winter_slope_median > 0 & low_25_winter_slope_lower > 0 ~ "Increase",
  	low_25_winter_slope_median < 0 & low_25_winter_slope_upper < 0 ~ "Decrease",
  	low_25_winter_slope_median > 0 & low_25_winter_slope_lower < 0 ~ "NotSig",
	low_25_winter_slope_median < 0 & low_25_winter_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_winter_50 = case_when(
  	med_50_winter_slope_median > 0 & med_50_winter_slope_lower > 0 ~ "Increase",
  	med_50_winter_slope_median < 0 & med_50_winter_slope_upper < 0 ~ "Decrease",
  	med_50_winter_slope_median > 0 & med_50_winter_slope_lower < 0 ~ "NotSig",
	med_50_winter_slope_median < 0 & med_50_winter_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_winter_IQR = case_when(
  	low_25_winter_slope_median < high_75_winter_slope_median ~ "Increase",
  	low_25_winter_slope_median > high_75_winter_slope_median ~ "Decrease",
  	.default = NA_character_
	)) %>%
	
	mutate(change_spring_75 = case_when(
  	high_75_spring_slope_median > 0 & high_75_spring_slope_lower > 0 ~ "Increase",
  	high_75_spring_slope_median < 0 & high_75_spring_slope_upper < 0 ~ "Decrease",
  	high_75_spring_slope_median > 0 & high_75_spring_slope_lower < 0 ~ "NotSig",
	high_75_spring_slope_median < 0 & high_75_spring_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_spring_25 = case_when(
  	low_25_spring_slope_median > 0 & low_25_spring_slope_lower > 0 ~ "Increase",
  	low_25_spring_slope_median < 0 & low_25_spring_slope_upper < 0 ~ "Decrease",
  	low_25_spring_slope_median > 0 & low_25_spring_slope_lower < 0 ~ "NotSig",
	low_25_spring_slope_median < 0 & low_25_spring_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_spring_50 = case_when(
  	med_50_spring_slope_median > 0 & med_50_spring_slope_lower > 0 ~ "Increase",
  	med_50_spring_slope_median < 0 & med_50_spring_slope_upper < 0 ~ "Decrease",
  	med_50_spring_slope_median > 0 & med_50_spring_slope_lower < 0 ~ "NotSig",
	med_50_spring_slope_median < 0 & med_50_spring_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_spring_IQR = case_when(
  	low_25_spring_slope_median < high_75_spring_slope_median ~ "Increase",
  	low_25_spring_slope_median > high_75_spring_slope_median ~ "Decrease",
  	.default = NA_character_
	)) %>%
	
	mutate(change_summer_75 = case_when(
  	high_75_summer_slope_median > 0 & high_75_summer_slope_lower > 0 ~ "Increase",
  	high_75_summer_slope_median < 0 & high_75_summer_slope_upper < 0 ~ "Decrease",
  	high_75_summer_slope_median > 0 & high_75_summer_slope_lower < 0 ~ "NotSig",
	high_75_summer_slope_median < 0 & high_75_summer_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_summer_25 = case_when(
  	low_25_summer_slope_median > 0 & low_25_summer_slope_lower > 0 ~ "Increase",
  	low_25_summer_slope_median < 0 & low_25_summer_slope_upper < 0 ~ "Decrease",
  	low_25_summer_slope_median > 0 & low_25_summer_slope_lower < 0 ~ "NotSig",
	low_25_summer_slope_median < 0 & low_25_summer_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_summer_50 = case_when(
  	med_50_summer_slope_median > 0 & med_50_summer_slope_lower > 0 ~ "Increase",
  	med_50_summer_slope_median < 0 & med_50_summer_slope_upper < 0 ~ "Decrease",
  	med_50_summer_slope_median > 0 & med_50_summer_slope_lower < 0 ~ "NotSig",
	med_50_summer_slope_median < 0 & med_50_summer_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_summer_IQR = case_when(
  	low_25_summer_slope_median < high_75_summer_slope_median ~ "Increase",
  	low_25_summer_slope_median > high_75_summer_slope_median ~ "Decrease",
  	.default = NA_character_
	)) %>%

	mutate(change_autumn_75 = case_when(
  	high_75_autumn_slope_median > 0 & high_75_autumn_slope_lower > 0 ~ "Increase",
  	high_75_autumn_slope_median < 0 & high_75_autumn_slope_upper < 0 ~ "Decrease",
  	high_75_autumn_slope_median > 0 & high_75_autumn_slope_lower < 0 ~ "NotSig",
	high_75_autumn_slope_median < 0 & high_75_autumn_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_autumn_25 = case_when(
  	low_25_autumn_slope_median > 0 & low_25_autumn_slope_lower > 0 ~ "Increase",
  	low_25_autumn_slope_median < 0 & low_25_autumn_slope_upper < 0 ~ "Decrease",
  	low_25_autumn_slope_median > 0 & low_25_autumn_slope_lower < 0 ~ "NotSig",
	low_25_autumn_slope_median < 0 & low_25_autumn_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_autumn_50 = case_when(
  	med_50_autumn_slope_median > 0 & med_50_autumn_slope_lower > 0 ~ "Increase",
  	med_50_autumn_slope_median < 0 & med_50_autumn_slope_upper < 0 ~ "Decrease",
  	med_50_autumn_slope_median > 0 & med_50_autumn_slope_lower < 0 ~ "NotSig",
	med_50_autumn_slope_median < 0 & med_50_autumn_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_autumn_IQR = case_when(
  	low_25_autumn_slope_median < high_75_autumn_slope_median ~ "Increase",
  	low_25_autumn_slope_median > high_75_autumn_slope_median ~ "Decrease",
  	.default = NA_character_
	)) %>%

	mutate(change_annual_75 = case_when(
  	high_75_annual_slope_median > 0 & high_75_annual_slope_lower > 0 ~ "Increase",
  	high_75_annual_slope_median < 0 & high_75_annual_slope_upper < 0 ~ "Decrease",
  	high_75_annual_slope_median > 0 & high_75_annual_slope_lower < 0 ~ "NotSig",
	high_75_annual_slope_median < 0 & high_75_annual_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_annual_25 = case_when(
  	low_25_annual_slope_median > 0 & low_25_annual_slope_lower > 0 ~ "Increase",
  	low_25_annual_slope_median < 0 & low_25_annual_slope_upper < 0 ~ "Decrease",
  	low_25_annual_slope_median > 0 & low_25_annual_slope_lower < 0 ~ "NotSig",
	low_25_annual_slope_median < 0 & low_25_annual_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_annual_50 = case_when(
  	med_50_annual_slope_median > 0 & med_50_annual_slope_lower > 0 ~ "Increase",
  	med_50_annual_slope_median < 0 & med_50_annual_slope_upper < 0 ~ "Decrease",
  	med_50_annual_slope_median > 0 & med_50_annual_slope_lower < 0 ~ "NotSig",
	med_50_annual_slope_median < 0 & med_50_annual_slope_upper > 0 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_annual_IQR = case_when(
  	low_25_annual_slope_median < high_75_annual_slope_median ~ "Increase",
  	low_25_annual_slope_median > high_75_annual_slope_median ~ "Decrease",
  	.default = NA_character_
	))



### Join the new_df back with mw_basins
new_df <- mw_basins %>% full_join(model_df, by = "STAID")

new_plot_df <- plot_df %>% full_join(mw_basins, by = "STAID") %>% drop_na(year)

### Write to file
write.csv(new_df, paste0(data_path, '/seasonal_runoff_trend_parameters.csv')) 
write.csv(new_plot_df, paste0(data_path, '/seasonal_runoff_trend_percentiles.csv'))


### Create an sf object of gauge locations
new_sf = st_as_sf(new_df, coords = c("LNG_GAGE", "LAT_GAGE"), crs = 4326 )
	+ scale_colour_manual(name = "Change", values = c("blue", "black",  "red"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
	+ theme_bw(12)
	
new_plot_sf = st_as_sf(new_plot_df, coords = c("LNG_GAGE", "LAT_GAGE"), crs = 4326 )
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


###########################################################################
###  Make plots for trends
###########################################################################

### WINTER: significant trends
p <- ggplot(data = new_sf %>% drop_na(change_winter_loc)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_winter_loc, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Winter Trend\nLocation", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_location_winter.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_location_winter.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_winter_scale)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_winter_scale, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Winter Trend\nScale", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_scale_winter.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_scale_winter.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_winter_shape)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_winter_shape, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Winter Trend\nShape", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_shape_winter.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_shape_winter.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_winter_freq)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_winter_freq, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Winter Trend\nFrequency", values = c("#19b0b3", "grey75",  "#b39d0e"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_frequency_winter.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_frequency_winter.svg"), p,  width = 6.5, height = 4.8)

max_value <- max(abs(new_sf$freq_winter), na.rm=TRUE)
p <- ggplot(data = new_sf %>% drop_na(freq_winter)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=freq_winter, shape=CLASS), size = 3) %>%
	+ scale_colour_distiller(name = "Winter Change\nFrequency\n(days/yr/yr)", palette = "RdBu", limits = c(-max_value, max_value), direction = 1, oob=squish) %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_frequency_winter_slope.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_frequency_winter_slope.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_winter_75, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Winter Trend\nInc Limb\nHigh 75%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_high75_winter.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_high75_winter.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_winter_25, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Winter Trend\nInc Limb\nLow 25%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_low25_winter.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_low25_winter.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_winter_50, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Winter Trend\nInc Limb\nMed 50%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_med50_winter.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_med50_winter.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_winter_IQR, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Winter Trend\nInc Limb\nIQR", values = c("#1f78b4", "#e31a1c"), breaks = c("Increase", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_IQR_winter.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_IQR_winter.svg"), p,  width = 6.5, height = 4.8)



### SPRING: significant trends
p <- ggplot(data = new_sf %>% drop_na(change_spring_loc)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_spring_loc, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Spring Trend\nLocation", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_location_spring.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_location_spring.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_spring_scale)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_spring_scale, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Spring Trend\nScale", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_scale_spring.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_scale_spring.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_spring_shape)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_spring_shape, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Spring Trend\nShape", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_shape_spring.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_shape_spring.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_spring_freq)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_spring_freq, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Spring Trend\nFrequency", values = c("#19b0b3", "grey75",  "#b39d0e"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_frequency_spring.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_frequency_spring.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_spring_75, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Spring Trend\nInc Limb\nHigh 75%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_high75_spring.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_high75_spring.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_spring_25, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Spring Trend\nInc Limb\nLow 25%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_low25_spring.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_low25_spring.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_spring_50, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Spring Trend\nInc Limb\nMed 50%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_med50_spring.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_med50_spring.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_spring_IQR, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Spring Trend\nInc Limb\nIQR", values = c("#1f78b4", "#e31a1c"), breaks = c("Increase", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_IQR_spring.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_IQR_spring.svg"), p,  width = 6.5, height = 4.8)



### SUMMER: significant trends
p <- ggplot(data = new_sf %>% drop_na(change_summer_loc)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_summer_loc, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Summer Trend\nLocation", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_location_summer.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_location_summer.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_summer_scale)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_summer_scale, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Summer Trend\nScale", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_scale_summer.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_scale_summer.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_summer_shape)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_summer_shape, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Summer Trend\nShape", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_shape_summer.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_shape_summer.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_summer_freq)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_summer_freq, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Summer Trend\nFrequency", values = c("#19b0b3", "grey75",  "#b39d0e"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_frequency_summer.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_frequency_summer.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_summer_75, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Summer Trend\nInc Limb\nHigh 75%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_high75_summer.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_high75_summer.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_summer_25, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Summer Trend\nInc Limb\nLow 25%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_low25_summer.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_low25_summer.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_summer_50, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Summer Trend\nInc Limb\nMed 50%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_med50_summer.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_med50_summer.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_summer_IQR, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Summer Trend\nInc Limb\nIQR", values = c("#1f78b4", "#e31a1c"), breaks = c("Increase", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_IQR_summer.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_IQR_summer.svg"), p,  width = 6.5, height = 4.8)




### AUTUMN: significant trends
p <- ggplot(data = new_sf %>% drop_na(change_autumn_loc)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_autumn_loc, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Autumn Trend\nLocation", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_location_autumn.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_location_autumn.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_autumn_scale)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_autumn_scale, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Autumn Trend\nScale", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_scale_autumn.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_scale_autumn.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_autumn_shape)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_autumn_shape, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Autumn Trend\nShape", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_shape_autumn.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_shape_autumn.svg"), p,  width = 6.5, height = 4.8)


p <- ggplot(data = new_sf %>% drop_na(change_autumn_freq)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_autumn_freq, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Autumn Trend\nFrequency", values = c("#19b0b3", "grey75",  "#b39d0e"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_frequency_autumn.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_frequency_autumn.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_autumn_75, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Autumn Trend\nInc Limb\nHigh 75%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_high75_autumn.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_high75_autumn.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_autumn_25, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Autumn Trend\nInc Limb\nLow 25%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_low25_autumn.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_low25_autumn.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_autumn_50, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Autumn Trend\nInc Limb\nMed 50%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_med50_autumn.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_med50_autumn.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_autumn_IQR, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Autumn Trend\nInc Limb\nIQR", values = c("#1f78b4", "#e31a1c"), breaks = c("Increase", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_IQR_autumn.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_IQR_autumn.svg"), p,  width = 6.5, height = 4.8)




### ANNUAL: significant trends
p <- ggplot(data = new_sf %>% drop_na(change_annual_loc)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_annual_loc, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Annual Trend\nLocation", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_location_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_location_annual.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_annual_scale)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_annual_scale, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Annual Trend\nScale", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_scale_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_scale_annual.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_annual_shape)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_annual_shape, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Annual Trend\nShape", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_shape_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_shape_annual.svg"), p,  width = 6.5, height = 4.8)


p <- ggplot(data = new_sf %>% drop_na(change_annual_freq)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_annual_freq, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Annual Trend\nFrequency", values = c("#19b0b3", "grey75",  "#b39d0e"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_frequency_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_frequency_annual.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_annual_75, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Annual Trend\nInc Limb\nHigh 75%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_high75_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_high75_annual.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_annual_25, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Annual Trend\nInc Limb\nLow 25%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_low25_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_low25_annual.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_annual_50, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Annual Trend\nInc Limb\nMed 50%", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_med50_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_med50_annual.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_plot_sf) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_annual_IQR, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Annual Trend\nInc Limb\nIQR", values = c("#1f78b4", "#e31a1c"), breaks = c("Increase", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_IQR_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_IQR_annual.svg"), p,  width = 6.5, height = 4.8)


