
# *------------------------------------------------------------------
# | FILE NAME: 01_streamflow_seasonal_gam.R
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
write_figures_path <- file.path(output_path, "figures_01_streamflow_seasonal_glm")
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



###########################################################################
###  GLM trend significance test - Mean, Shape
###########################################################################
parameterCd <- "00060"  # Discharge

### Create empty data frame
id_list <- mw_basins$STAID

model_df <- data.frame(STAID = id_list) %>%
	mutate(mean_winter = NA) %>%
	mutate(mean_spring = NA) %>%
	mutate(mean_summer = NA) %>%
	mutate(mean_autumn = NA) %>%
	mutate(mean_annual = NA) %>%
	mutate(shape_winter = NA) %>%
	mutate(shape_spring = NA) %>%
	mutate(shape_summer = NA) %>%
	mutate(shape_autumn = NA) %>%
	mutate(shape_annual = NA) %>%
	mutate(p_mean_winter = NA) %>%
	mutate(p_mean_spring = NA) %>%
	mutate(p_mean_summer = NA) %>%
	mutate(p_mean_autumn = NA) %>%
	mutate(p_mean_annual = NA) %>%
	mutate(p_shape_winter = NA) %>%
	mutate(p_shape_spring = NA) %>%
	mutate(p_shape_summer = NA) %>%
	mutate(p_shape_autumn = NA) %>%
	mutate(p_shape_annual = NA)
	


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
	  
	### Determine drainage area (sq km)
	drain_area_sqkm <- mw_basins$DRAIN_SQKM[k]
	
	### Convert units of flow (ft3/s) to runoff (mm/day)
	flow_k <- flow_k %>% 
		rename(flow_cfs = flow) %>%	
		mutate(flow_mm_day = flow_cfs / drain_area_sqkm / (3280.84)^3 * 10^6 * 86400)
	  
	#ggplot(flow_k, aes(x=year, y=flow_mm_day)) + geom_point() + scale_y_continuous(name = "Daily Discharge (mm per unit area)") + theme_classic(12)
	#ggplot(flow_k, aes(x=flow_cfs)) + geom_histogram() + theme_classic(12)
	#ggplot(flow_k, aes(x=year, y=flow_cfs)) + geom_point()+ theme_classic(12)



	### Filter for seasons
	flow_k_winter <- flow_k%>%
		filter(month == "1" | month == "2" | month == "12") %>%
		as.data.frame() %>% 
		mutate(year = case_when(month == 12 ~ year+1, 
					   TRUE  ~ year))

	flow_k_spring <- flow_k%>%
		filter(month == "3" | month == "4" | month == "5") %>%
		as.data.frame()

	flow_k_summer <- flow_k%>%
		filter(month == "6" | month == "7" | month == "8") %>%
		as.data.frame()

	flow_k_autumn <- flow_k%>%
		filter(month == "9" | month == "10" | month == "11") %>%
		as.data.frame()
	
	flow_k_annual <- flow_k

	### Calculate seasonal mean/median
	flow_k_winter_avg <- flow_k_winter %>%
		group_by(year) %>%
		summarize(mean_flow = mean(flow_mm_day, na.rm=TRUE), median_flow=median(flow_mm_day), count = n())

	flow_k_spring_avg <- flow_k_spring %>%
		group_by(year) %>%
		summarize(mean_flow = mean(flow_mm_day, na.rm=TRUE), median_flow=median(flow_mm_day), count = n())		

	flow_k_summer_avg <- flow_k_summer %>%
		group_by(year) %>%
		summarize(mean_flow = mean(flow_mm_day, na.rm=TRUE), median_flow=median(flow_mm_day), count = n())		

	flow_k_autumn_avg <- flow_k_autumn %>%
		group_by(year) %>%
		summarize(mean_flow = mean(flow_mm_day, na.rm=TRUE), median_flow=median(flow_mm_day), count = n())
	
	flow_k_annual_avg <- flow_k_annual %>%
		group_by(year) %>%
		summarize(mean_flow = mean(flow_mm_day, na.rm=TRUE), median_flow=median(flow_mm_day), count = n())


	### Filter to only years with more than 75 days of observations per season (no more than 2 weeks missing)
	flow_k_winter_avg <- flow_k_winter_avg %>% filter(count > 75)
	flow_k_spring_avg <- flow_k_spring_avg %>% filter(count > 75)
	flow_k_summer_avg <- flow_k_summer_avg %>% filter(count > 75)
	flow_k_autumn_avg <- flow_k_autumn_avg %>% filter(count > 75)


	### Make sure there are at least 70 years of flow data for each season
	if(dim(flow_k_winter_avg)[1] > 70 & dim(flow_k_spring_avg)[1] > 70 & dim(flow_k_summer_avg)[1] > 70 & dim(flow_k_autumn_avg)[1] > 70){


		### Make sure median seasonal flows > 0 (no intermittent or ephemeral streams)
		if(((min(flow_k_winter_avg$median_flow)) * (min(flow_k_spring_avg$median_flow)) * (min(flow_k_summer_avg$median_flow)) * (min(flow_k_autumn_avg$median_flow))) > 0) {
		
		### Make sure there is recent flow data (at least through 2010)
		if(max(flow_k_winter_avg$year) > 2009){

		###########################################################################
		###  Fit a GAM for each season
		###########################################################################
		### Fit gammals model
		gammals_fit_winter <- gam(list(mean_flow ~ year,
			~ year),
			family=gammals(link=list("identity","identity")), data = flow_k_winter_avg,
			select = TRUE)

		gammals_fit_spring <- gam(list(mean_flow ~ year,
			~ year),
			family=gammals(link=list("identity","identity")), data = flow_k_spring_avg,
			select = TRUE)

		gammals_fit_summer <- gam(list(mean_flow ~ year,
			~ year),
			family=gammals(link=list("identity","identity")), data = flow_k_summer_avg,
			select = TRUE)

		gammals_fit_autumn <- gam(list(mean_flow ~ year,
			~ year),
			family=gammals(link=list("identity","identity")), data = flow_k_autumn_avg,
			select = TRUE)

		gammals_fit_annual <- gam(list(mean_flow ~ year,
			~ year),
			family=gammals(link=list("identity","identity")), data = flow_k_annual_avg,
			select = TRUE)



		### Extract results
		mean_winter <- summary(gammals_fit_winter)$p.coeff["year"]
		shape_winter <- summary(gammals_fit_winter)$p.coeff["year.1"]
		p_mean_winter <- summary(gammals_fit_winter)$p.pv["year"]
		p_shape_winter <- summary(gammals_fit_winter)$p.pv["year.1"]
		
		mean_spring <- summary(gammals_fit_spring)$p.coeff["year"]
		shape_spring <- summary(gammals_fit_spring)$p.coeff["year.1"]
		p_mean_spring <- summary(gammals_fit_spring)$p.pv["year"]
		p_shape_spring <- summary(gammals_fit_spring)$p.pv["year.1"]

		mean_summer <- summary(gammals_fit_summer)$p.coeff["year"]
		shape_summer <- summary(gammals_fit_summer)$p.coeff["year.1"]
		p_mean_summer <- summary(gammals_fit_summer)$p.pv["year"]
		p_shape_summer <- summary(gammals_fit_summer)$p.pv["year.1"]
		
		mean_autumn <- summary(gammals_fit_autumn)$p.coeff["year"]
		shape_autumn <- summary(gammals_fit_autumn)$p.coeff["year.1"]
		p_mean_autumn <- summary(gammals_fit_autumn)$p.pv["year"]
		p_shape_autumn <- summary(gammals_fit_autumn)$p.pv["year.1"]
		
		mean_annual <- summary(gammals_fit_annual)$p.coeff["year"]
		shape_annual <- summary(gammals_fit_annual)$p.coeff["year.1"]
		p_mean_annual <- summary(gammals_fit_annual)$p.pv["year"]
		p_shape_annual <- summary(gammals_fit_annual)$p.pv["year.1"]


		
		### Save into columns
		model_df$mean_winter[k] <- mean_winter
		model_df$shape_winter[k] <- shape_winter
		model_df$p_mean_winter[k] <- p_mean_winter
		model_df$p_shape_winter[k] <- p_shape_winter
		
		model_df$mean_spring[k] <- mean_spring
		model_df$shape_spring[k] <- shape_spring
		model_df$p_mean_spring[k] <- p_mean_spring
		model_df$p_shape_spring[k] <- p_shape_spring

		model_df$mean_summer[k] <- mean_summer
		model_df$shape_summer[k] <- shape_summer
		model_df$p_mean_summer[k] <- p_mean_summer
		model_df$p_shape_summer[k] <- p_shape_summer

		model_df$mean_autumn[k] <- mean_autumn
		model_df$shape_autumn[k] <- shape_autumn
		model_df$p_mean_autumn[k] <- p_mean_autumn
		model_df$p_shape_autumn[k] <- p_shape_autumn

		model_df$mean_annual[k] <- mean_annual
		model_df$shape_annual[k] <- shape_annual
		model_df$p_mean_annual[k] <- p_mean_annual
		model_df$p_shape_annual[k] <- p_shape_annual



	}
	}
	}

}, error = function(e) {
	### Log error and continue
	cat("Error encountered for gauge", k, ": ", e$message, "\n")
	})

}


### Add columns that do if then statements based on change and p-value
model_df <- model_df %>%
	mutate(change_winter_mean = NA) %>%
	mutate(change_spring_mean = NA) %>%
	mutate(change_summer_mean = NA) %>%
	mutate(change_autumn_mean = NA) %>%
	mutate(change_annual_mean = NA) %>%
	mutate(change_winter_mean = case_when(
  	mean_winter > 0 & p_mean_winter < 0.05 ~ "Increase",
  	mean_winter < 0 & p_mean_winter < 0.05 ~ "Decrease",
  	p_mean_winter >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_spring_mean = case_when(
  	mean_spring > 0 & p_mean_spring < 0.05 ~ "Increase",
  	mean_spring < 0 & p_mean_spring < 0.05 ~ "Decrease",
  	p_mean_spring >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_summer_mean = case_when(
  	mean_summer > 0 & p_mean_summer < 0.05 ~ "Increase",
  	mean_summer < 0 & p_mean_summer < 0.05 ~ "Decrease",
  	p_mean_summer >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_autumn_mean = case_when(
  	mean_autumn > 0 & p_mean_autumn < 0.05 ~ "Increase",
  	mean_autumn < 0 & p_mean_autumn < 0.05 ~ "Decrease",
  	p_mean_autumn >= 0.05 ~ "NotSig",
  	.default = NA_character_
	)) %>%
	mutate(change_annual_mean = case_when(
  	mean_annual > 0 & p_mean_annual < 0.05 ~ "Increase",
  	mean_annual < 0 & p_mean_annual < 0.05 ~ "Decrease",
  	p_mean_annual >= 0.05 ~ "NotSig",
  	.default = NA_character_
	))



model_df <- model_df %>%
	mutate(change_winter_shape = NA) %>%
	mutate(change_spring_shape = NA) %>%
	mutate(change_summer_shape = NA) %>%
	mutate(change_autumn_shape = NA) %>%
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




### Join the new_df back with mw_basins
new_df <- mw_basins %>% full_join(model_df, by = "STAID")
new_df <- new_df %>% select(-'HCDN-2009') %>% drop_na(change_winter_mean)

### Join with filtered storage df
low_storage_df <- read_csv(paste0(data_path, '/upstream_storage_low.csv'))

new_df <- new_df %>%
	filter(STAID %in% low_storage_df$STAID) %>%
	filter(STAID != "05444000") %>% filter(STAID != "04094000") %>% filter(STAID != "03301500")


### Write to file. 
write.csv(new_df, paste0(data_path, '/seasonal_flow_trends.csv'))


### Create an sf object of gauge locations
new_sf = st_as_sf(new_df, coords = c("LNG_GAGE", "LAT_GAGE"), crs = 4326 )
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

### Winter: significant trends
p <- ggplot(data = new_sf %>% drop_na(change_winter_mean)) %>%
#p <- ggplot(data = new_sf %>% drop_na(change_winter_mean) %>% filter(CLASS=='Ref')) %>%
		#+ geom_stars(data = background, alpha = 0.7) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_winter_mean, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Winter Trend\nMean", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_mean_winter.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_mean_winter.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_winter_shape)) %>%
		#+ geom_stars(data = background, alpha = 0.7) %>%
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


### Spring: significant trends
p <- ggplot(data = new_sf %>% drop_na(change_spring_mean)) %>%
#p <- ggplot(data = new_sf %>% drop_na(change_spring_mean) %>% filter(CLASS=='Ref')) %>%
		#+ geom_stars(data = background, alpha = 0.7) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_spring_mean, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Spring Trend\nMean", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_mean_spring.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_mean_spring.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_spring_shape)) %>%
		#+ geom_stars(data = background, alpha = 0.7) %>%
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


### Summer: significant trends
p <- ggplot(data = new_sf %>% drop_na(change_summer_mean)) %>%
#p <- ggplot(data = new_sf %>% drop_na(change_summer_mean) %>% filter(CLASS=='Ref')) %>%
		#+ geom_stars(data = background, alpha = 0.7) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_summer_mean, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Summer Trend\nMean", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_mean_summer.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_mean_summer.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_summer_shape)) %>%
		#+ geom_stars(data = background, alpha = 0.7) %>%
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


### Autumn: significant trends
p <- ggplot(data = new_sf %>% drop_na(change_autumn_mean)) %>%
#p <- ggplot(data = new_sf %>% drop_na(change_autumn_mean) %>% filter(CLASS=='Ref')) %>%
		#+ geom_stars(data = background, alpha = 0.7) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_autumn_mean, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Autumn Trend\nMean", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_mean_autumn.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_mean_autumn.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_autumn_shape)) %>%
		#+ geom_stars(data = background, alpha = 0.7) %>%
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


### Annual: significant trends
p <- ggplot(data = new_sf %>% drop_na(change_annual_mean)) %>%
#p <- ggplot(data = new_sf %>% drop_na(change_annual_mean) %>% filter(CLASS=='Ref')) %>%
		#+ geom_stars(data = background, alpha = 0.7) %>%
		+ geom_sf(data = lakes, fill = "#9bbff4") %>%
		+ geom_sf(data = rivers, colour = "#4a80f5") %>%
		+ geom_sf(data = usa_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(data = ca_states, fill = NA, alpha = 0.5) %>%
		+ geom_sf(aes(color=change_annual_mean, shape=CLASS), size = 3) %>%
		+ scale_colour_manual(name = "Annual Trend\nMean", values = c("#1f78b4", "grey75",  "#e31a1c"), breaks = c("Increase",  "NotSig", "Decrease"), na.value="grey80") %>%
		+ scale_fill_identity() %>%
		+ theme_bw(12) %>%
		+ theme(legend.position = "right") %>%
		+ theme( panel.grid.minor = element_blank(), panel.grid.major = element_blank()) %>%
		+ coord_sf(xlim = c(-93, -77.5), ylim = c(36.4, 47.5), expand = FALSE) %>%
		+ xlab("Longitude") %>%
		+ ylab("Latitude")
	### Save plot
	ggsave(file.path(write_figures_path, "map_trend_mean_annual.png"), p,  width = 6.5, height = 4.8, dpi = 600)
	ggsave(file.path(write_figures_path, "map_trend_mean_annual.svg"), p,  width = 6.5, height = 4.8)

p <- ggplot(data = new_sf %>% drop_na(change_annual_shape)) %>%
		#+ geom_stars(data = background, alpha = 0.7) %>%
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




