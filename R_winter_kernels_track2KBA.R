#### --------------------------------------------------------------------------~
#### WINTER DISTRIBUTION OF ADÉLIE PENGUINS                         ------------
#### KERNEL UTILIZATION DISTRIBUTION
#### REPRESENTATIVNESS {track2KBA}
#### SITE OF POTENTIAL IMPORTANCE FOR CONSERVATION {track2KBA}
#### --------------------------------------------------------------------------~

#### Code for manuscript:
#### "Do penguins care about on-paper boundaries?
#### Conservation implications of spatio-temporal winter consistency 
#### in an Antarctic sentinel species"

#### Zuzana Zajková
#### zajkova@icm.csic.es

#### --------------------------------------------------------------------------~
#### LIBRARIES ---- 
#### --------------------------------------------------------------------------~

my_path <- "your_path"

library(tidyverse)
library(sf)
library(sp)
library(track2KBA) ## version 1.1.2
library(here)

library(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf")


#### --------------------------------------------------------------------------~
#### 1) GETTING DATA ----
#### --------------------------------------------------------------------------~

#### Data previously downloaded from Movebank and transformed to sf format
#### https://www.movebank.org/, study ID: 4778411162

load("Adelies_GLS_DDU_movebank.RData")

dim(tracks_sf)
# 23986     9


tracks_sf <- 
  tracks_sf %>% 
  dplyr::filter(stage == "wintering") %>% 
  dplyr::mutate(lon = as.numeric(longitude),
                lat = as.numeric(latitude), 
                month = lubridate::month(timestamp, label = FALSE), 
                month_label = lubridate::month(timestamp, label = TRUE, abbr = FALSE), 
                winter_year = as.character(lubridate::year(timestamp)), 
                gls_id = track_id,
                pit_id = bird_id,
                winter_year_gls_id = paste(winter_year, gls_id, sep ="_")) %>% 
  dplyr::select(winter_year, winter_year_gls_id, gls_id, 
                lon, lat, 
                date_gmt, time_gmt, month, month_label)


tracks_df <- 
  tracks_sf %>% 
  as.data.frame()


dim(tracks_sf)
# 21733     9
dim(tracks_df)

length(unique(tracks_df$winter_year_gls_id))
# 61 (excluded winter_year_gls_id 3031)

#### --------------------------------------------------------------------------~
#### --------------------------------------------------------------------------~
#### 2) KERNEL UD ----
#### --------------------------------------------------------------------------~
#### --------------------------------------------------------------------------~

#### USING ALL 61 TRACKS, INLCUDING NOT-COMPETE TRIPS
sort(table(tracks_df$winter_year_gls_id))

#### Convert to Lambert Azimuthal Equal Area for kernel calculations
mid_point <- data.frame(geosphere::centroid(cbind(tracks_sp$lon, 
                                                  tracks_sp$lat)))
mid_point
# lon       lat
# 123.0169 -63.77194

laea_proj <- paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep="")

tracks_sp_laea <- spTransform(tracks_sp, 
                              CRS = laea_proj) 

plot(tracks_sp_laea)

#### --------------------------------------------------------------------------~
#### KERNEL CALCULATIONS PER YEAR ----
#### --------------------------------------------------------------------------~

library(adehabitatHR)

#### Calculate kernels by year (year needs to be the first column)
KUD_YEAR <- kernelUD(tracks_sp_laea[, 1], h = 186000, same4all = TRUE, grid = 1000)

#### Save 
save(KUD_YEAR, 
     file = "KUD_YEAR_final.RData")

image(KUD_YEAR)

#### Get vertices
KUD_YEAR5  <- getverticeshr(KUD_YEAR,  5)
KUD_YEAR25 <- getverticeshr(KUD_YEAR, 25)
KUD_YEAR50 <- getverticeshr(KUD_YEAR, 50)
KUD_YEAR75 <- getverticeshr(KUD_YEAR, 75)
KUD_YEAR95 <- getverticeshr(KUD_YEAR, 95)

#### Exploration plots
plot(KUD_YEAR95)
plot(KUD_YEAR75, col = 2, add = TRUE)
plot(KUD_YEAR50, col = 3, add = TRUE)
plot(KUD_YEAR25, col = 4, add = TRUE)
plot(KUD_YEAR5,  col = 5, add = TRUE)

KUD_YEAR50_sf <-
  st_as_sf(KUD_YEAR50) %>% 
  rename(winter_year = id) %>% 
  mutate(kud = "kud50")

KUD_YEAR95_sf <- 
  st_as_sf(KUD_YEAR95) %>% 
  rename(winter_year = id)%>% 
  mutate(kud = "kud95")


#### Plot 
ggplot() + 
  geom_sf(data = world) +
  ## KUD
  geom_sf(data = KUD_YEAR95_sf, aes(fill = winter_year), colour = NA, alpha = 0.1) +
  geom_sf(data = KUD_YEAR50_sf, aes(fill = winter_year, colour = winter_year), alpha = 0.2) +
  scale_colour_viridis_d(option = "magma", begin = 0.25, name = "Year") +
  scale_fill_viridis_d(option = "magma", begin = 0.25, name = "Year") +
  coord_sf() 


#### --------------------------------------------------------------------------~
#### --------------------------------------------------------------------------~
#### 3) track2KBA ----
#### --------------------------------------------------------------------------~
#### --------------------------------------------------------------------------~

# Beal, M., Oppel, S., Handley, J., Pearmain, E. J., Morera‐Pujol, V., 
# Carneiro, A. P., ... & Dias, M. P. (2021). 
# track2KBA: An R package for identifying important sites for biodiversity 
# from tracking data. Methods in Ecology and Evolution.
# https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.13713


#### --------------------------------------------------------------------------~
#### REPRESENTATIVNESS ---- RUN LOOP OVER YEARS  ----
#### --------------------------------------------------------------------------~


#### FOR 50% KDE ---
years <- c(unique(tracks_df$winter_year), "all")

repr_years_50 <- data.frame()


set.seed(22)

for(i in seq_along(years)) {
  
  sel_year <- years[i]
  
  if (sel_year != "all") {
    tracks_df_sel <- filter(tracks_df, winter_year == sel_year)
    message(paste("Calculating representativeness for year", sel_year), ".")
  }
  
  if (sel_year == "all") {
    tracks_df_sel <- tracks_df
    message(paste("Calculating representativeness for all years pooled."))
  }
  
  dataGroup <- formatFields(dataGroup = tracks_df_sel, 
                            fieldID   = "winter_year_gls_id", 
                            fieldDate = "date_gmt", 
                            fieldTime = "time_gmt",
                            fieldLon  = "lon", 
                            fieldLat  = "lat")
  
  # Split to trips
  trips <- tripSplit(dataGroup  = dataGroup,
                     colony     = colony,
                     innerBuff  = 100,             # kilometers
                     returnBuff = 100,
                     duration   = 30*24,           # hours
                     rmNonTrip  = FALSE)
  
  # # Summary of trip stats
  # sumTrips <- tripSummary(trips = trips, colony = colony)
  # sumTrips
  
  #### Transform the tracking data to an equal-area projection. 
  tracks <- projectTracks(dataGroup = trips, projType = 'azim', custom = TRUE)
  
  
  #### Kernel UD
  KDE_50 <- estSpaceUse(tracks = tracks, 
                        scale   = 186, ## assuming the spatial error in GLS (Phillips et al. 2004)
                        levelUD = 50, 
                        polyOut = TRUE)
  #### Plot
  mapKDE(KDE = KDE_50$UDPolygons, colony = colony) 
  
  ggsave(plot = last_plot(), 
         filename = paste0(my_path, "/PLOTS/OVERLAPS/PLOT_Representativeness_50KDE_tracks_final_", 
                           sel_year, ".png"), 
         width = 10, height = 6, units = "in", dpi = 300)
  
  #### Representativness of the sample
  #### Takes a while to calculate
  
  png(filename = paste0(my_path, "/PLOTS/OVERLAPS/PLOT_Representativeness_50KDE_tracks_final_asymptote_", 
                        sel_year, ".png"), 
      width = 10, height = 8, units = "in", res = 600)
  
  repr <- repAssess(tracks    = tracks, 
                    KDE       = KDE_50$KDE.Surface,
                    levelUD   = 50,
                    iteration = 100, 
                    bootTable = TRUE)
  
  dev.off()
  
  #### Save output
  repr_output <- repr[[1]]
  repr_output$winter_year <- sel_year
  repr_output$iterations  <- 100
  repr_years_50 <- dplyr::bind_rows(repr_years_50, repr_output)
  
  #### Save output table
  save(repr, file = paste0(my_path, "/RESULTS/ADELIES_repr50KDE_fulloutputtable_final_", sel_year, ".RData"))
  
}


#### Save output
save(repr_years_50, file = paste0(my_path, "/RESULTS/ADELIES_repr50KDE_years_final.RData"))


#### --------------------------------------------------------------------------~
#### SITE OF POTENTIAL IMPORTANCE FOR CONSERVATION ----
#### --------------------------------------------------------------------------~

### For individual kernels for 50% UD.
file = paste0(my_path, "/RESULTS/ADELIES_repr50KDE_years_final.RData")
repr <- 
  repr_years_50 %>%
  filter(winter_year == "all")

set.seed(22)

#### Density 
Site_density <- findSite(KDE = KDE_50$KDE.Surface, # calculated for all
                         represent = repr$out,
                         levelUD = 50,
                         popSize = 40000,  # N individual seabirds breed one the site (Barbraud et al., 2020)
                         polyOut = FALSE)

mapSite(Site_density, colony = colony) 

Site <- findSite(KDE = KDE_50$KDE.Surface,
                 represent = repr$out,
                 levelUD = 50,
                 popSize = 40000,     # N individual seabirds breed one the site (Barbraud et al., 2020)
                 polyOut = TRUE)

mapSite(Site, colony = colony) 

#### Save output
save(Site, file = paste0(my_path, "/RESULTS/ADELIES_repr50KDE_Site_final.RData"))


