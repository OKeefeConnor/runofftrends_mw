
# *------------------------------------------------------------------
# | FILE NAME: 00_storage.R
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
write_figures_path <- file.path(output_path, "figures_00_storage")
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
	filter(STATE == "OH" | STATE == "IN" | STATE == "KY" | STATE == "MI" | STATE == "WI" | STATE == "WV" | STATE == "IL") %>%
	as.data.frame()
	
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
		
		### Calculate daily mean and median flow
		flow_k_daily <- flow_k %>%
			summarize(mean_flow = mean(flow, na.rm=TRUE), median_flow=median(flow))	
		
		flow_k_annual_df <- as.data.frame(flow_k_annual) %>%
			mutate(STAID = gauge_k)

		flow_k_daily_df <- as.data.frame(flow_k_daily) %>%
			mutate(STAID = gauge_k)



	if (k == 1){
		flow_annual_df <- flow_k_annual_df
		
		flow_df <- flow_k_daily_df
		
	} else {
		flow_annual_df <- flow_annual_df %>%
			bind_rows(flow_k_annual_df)
		
		flow_df <- flow_df %>%
			bind_rows(flow_k_daily_df)		
	}


	
	#}
}, error = function(e) {
	### Log error and continue
	cat("Error encountered for gauge", k, ": ", e$message, "\n")
	})



}



### Combine data frames
flow_df <- flow_df %>% 
	full_join(mw_basins, by = "STAID") %>% 
	drop_na(mean_flow)

flow_annual_df <- flow_annual_df %>% 
	full_join(mw_basins, by = "STAID") %>% 
	drop_na(mean_flow_annual)



### Convert mean daily flow (cfs) to mean daily runoff (mm/day)
summary_df <- flow_df %>% 
	drop_na(mean_flow) %>%
	mutate(mean_flow_mm_day = mean_flow / DRAIN_SQKM / (3280.84)^3 * 10^6 * 86400)



#######################################
### Read in dam storage information ###
storage_df_1 <- read_csv(file.path(gages_folder, "Dataset1_BasinID/conterm_hydromod_dams.txt"))
hist(storage_df_1$STOR_NID_2009)

### Combine storage with gauge information
storage_df <- left_join(summary_df, storage_df_1, by = "STAID")

### Calculate number of days of storage for each gauge
storage_df <- storage_df %>%
		### megaliters per sq km to cubic feet
	mutate(storage_ft3 = (STOR_NID_2009 * DRAIN_SQKM * 1000 * 35.3147)) %>%
		### cfs to cubic feet per day
	mutate(flow_ft3_day = (mean_flow * 86400)) %>%
		### days of storage
	mutate(storage_days = (storage_ft3 / flow_ft3_day))

### Write to file. 
write.csv(storage_df, paste0(data_path, '/upstream_storage.csv'))



### Filter for gages with less than 48 days of storage
low_storage_df <- storage_df %>%
	filter(storage_days < 48)

### Write to file. 
write.csv(low_storage_df, paste0(data_path, '/upstream_storage_low.csv'))



### Create an sf object of gauge locations
plot_storage_sf = st_as_sf(storage_df, coords = c("LNG_GAGE", "LAT_GAGE"), crs = 4326 ) 
#	+ scale_colour_manual(name = "Change", values = c("blue", "black",  "red"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
#	+ theme_bw(12)

### Map data comes from https://www.naturalearthdata.com/downloads/
### Download state borders
usa_states <- ne_states(country = 'United States of America', returnclass = 'sf')
ca_states <- ne_states(country = 'Canada', returnclass = 'sf')
mx_states <- ne_states(country = 'Mexico', returnclass = 'sf')
### Download lakes and rivers
lakes <- ne_download(scale = 50, type = 'lakes', category = 'physical', returnclass = 'sf')
rivers <- ne_download(scale = 10, type = 'rivers_lake_centerlines', category = 'physical', returnclass = 'sf')


### Plot map of storage for gauges
p <- ggplot(data = plot_storage_sf %>% drop_na(storage_days)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=storage_days, shape=CLASS), size = 3) %>%
		#+ scale_colour_manual(name = "Average Flow\n(mm/day)", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_colour_distiller(name = "Days of\nStorage", palette = "Spectral", limits = c(0, 50), oob=squish, direction=-1) %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_storage.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_storage.svg"), p,  width = 6.5, height = 4.8)





################################################
### Read in basin classification information ###
basin_classif_1 <- read_csv(file.path(gages_folder, "Dataset1_BasinID/conterm_bas_classif.txt"))
basin_classif_1 <- as.data.frame(basin_classif_1) %>%
	select(-CLASS)
hist(basin_classif_1$HYDRO_DISTURB_INDX)

### Join with mw_basins
basin_classif <- left_join(flow_df, basin_classif_1, by = "STAID") %>%
	drop_na(HYDRO_DISTURB_INDX)
hist(basin_classif$HYDRO_DISTURB_INDX)



### Create an sf object of gauge locations
plot_sf = st_as_sf(basin_classif, coords = c("LNG_GAGE", "LAT_GAGE"), crs = 4326 )
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


### Plot map of hydrologic disturbance for gauges
p <- ggplot(data = plot_sf %>% drop_na(HYDRO_DISTURB_INDX)) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=HYDRO_DISTURB_INDX, shape=CLASS), size = 3) %>%
		#+ scale_colour_manual(name = "Average Flow\n(mm/day)", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_colour_distiller(name = "Hydrologic\nDisturbance\nIndex", palette = "Spectral", direction=-1) %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_hydro_disturb_indx.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_hydro_disturb_indx.svg"), p,  width = 6.5, height = 4.8)


