
# *------------------------------------------------------------------
# | FILE NAME: 05_gauge_trend_info.R
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
write_figures_path <- file.path(output_path, "figures_05_gauge_trend_info")
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

parameterCd <- "00060"  # Discharge



### Create empty data frame
id_list <- mw_basins$STAID

model_df <- data.frame(STAID = id_list)



### Run loop for each gauge
for (k in seq(1, dim(model_df)[1])){

cat(k)
cat("\n")

tryCatch({

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
	flow_k <- flow_k %>%
		add_column(year = year(flow_k$date))
	flow_k <- flow_k %>%
		add_column(month = month(flow_k$date))
	flow_k <- flow_k %>%
		add_column(jdate = day(flow_k$date))
	
	

		### Calculate annual mean and median
		flow_k_annual <- flow_k %>%
			group_by(year) %>%
			summarize(mean_flow_annual = mean(flow, na.rm=TRUE), median_flow_annual = median(flow, na.rm = TRUE))	
		
		flow_k_annual_df <- as.data.frame(flow_k_annual) %>%
			mutate(STAID = gauge_k)




	if (k == 1){
		flow_annual_df <- flow_k_annual_df
		
	} else {
		flow_annual_df <- flow_annual_df %>%
			bind_rows(flow_k_annual_df)

	}


	
}, error = function(e) {
	### Log error and continue
	cat("Error encountered for gauge", k, ": ", e$message, "\n")
	})



}



### Combine data frames
flow_annual_df <- flow_annual_df %>% 
	full_join(mw_basins, by = "STAID") %>% 
	drop_na(mean_flow_annual)



### Count the number of gages for each year
count_df <- flow_annual_df %>%
	group_by(year) %>%
	summarize(num_gages = n())

### Plot number of gauges per year
p <- ggplot(count_df, aes(x = year, y = num_gages)) %>%
	+ geom_bar(stat = "identity", fill = "black", alpha = 1) %>%
	+ labs(x = "Year", y = "Number of Gages")
	### Save plot
	ggsave(file.path(write_figures_path, "number_gages_year.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "number_gages_year.svg"), p,  width = 6.5, height = 4.8)







####################################################################
### Read in saved data frames to count # of inc/dec gauges
####################################################################

### Seasonal Streamflow: Mean
seasonal_flow_mean <- read_csv(paste0(data_path, '/seasonal_flow_trends.csv'))

sum(seasonal_flow_mean$change_winter_mean == 'Increase')
sum(seasonal_flow_mean$change_winter_mean == 'Decrease')
sum(seasonal_flow_mean$change_winter_mean == 'NotSig')

sum(seasonal_flow_mean$change_spring_mean == 'Increase')
sum(seasonal_flow_mean$change_spring_mean == 'Decrease')
sum(seasonal_flow_mean$change_spring_mean == 'NotSig')

sum(seasonal_flow_mean$change_summer_mean == 'Increase')
sum(seasonal_flow_mean$change_summer_mean == 'Decrease')
sum(seasonal_flow_mean$change_summer_mean == 'NotSig')

sum(seasonal_flow_mean$change_autumn_mean == 'Increase')
sum(seasonal_flow_mean$change_autumn_mean == 'Decrease')
sum(seasonal_flow_mean$change_autumn_mean == 'NotSig')



### Seasonal Runoff: Median
seasonal_runoff_median <- read_csv(paste0(data_path, '/seasonal_runoff_trend_percentiles.csv'))

seasonal_runoff_median <- seasonal_runoff_median %>%
  group_by(STANAME) %>%
  slice(1) %>%
  ungroup()

sum(seasonal_runoff_median$change_winter_50 == 'Increase')
sum(seasonal_runoff_median$change_winter_50 == 'Decrease')
sum(seasonal_runoff_median$change_winter_50 == 'NotSig')

sum(seasonal_runoff_median$change_spring_50 == 'Increase')
sum(seasonal_runoff_median$change_spring_50 == 'Decrease')
sum(seasonal_runoff_median$change_spring_50 == 'NotSig')

sum(seasonal_runoff_median$change_summer_50 == 'Increase')
sum(seasonal_runoff_median$change_summer_50 == 'Decrease')
sum(seasonal_runoff_median$change_summer_50 == 'NotSig')

sum(seasonal_runoff_median$change_autumn_50 == 'Increase')
sum(seasonal_runoff_median$change_autumn_50 == 'Decrease')
sum(seasonal_runoff_median$change_autumn_50 == 'NotSig')

