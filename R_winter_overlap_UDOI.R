#### --------------------------------------------------------------------------~
#### WINTER DISTRIBUTION OF ADÉLIE PENGUINS                         ------------
#### OVERLAP KERNEL UD
#### Utilization distribution overlap index (UDOI) 
#### RANDOMIZATION FOR STATISTICAL TESTING
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

library(tidyverse)
library(adehabitatHR)
library(sf)
library(sp)
# library(geosphere)

my_path <- "your_path"

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
#### 2) PROJECT DATA                                                        ----
#### --------------------------------------------------------------------------~

#### Lambert Azimuthal Equal Area
mid_point <- data.frame(geosphere::centroid(cbind(tracks_sp$lon, 
                                                  tracks_sp$lat)))
# lon       lat
# 123.0169 -63.77194

laea_proj <- paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep="")

tracks_sp_laea <- spTransform(tracks_sp, 
                              CRS = laea_proj) 

plot(tracks_sp_laea)

# best_tracks_sp_winter_filtered_laea 


#### --------------------------------------------------------------------------~
#### 3) RANDOMIZE THE YEAR AND CALCULATE THE OVERLAPS                       ----
#### --------------------------------------------------------------------------~

#### --------------------------------------------------------------------------~
#### CALCULATE THE RANDOMIZED VALUES OF OVERLAPS                            ---- 
#### --------------------------------------------------------------------------~

#### "we performed a randomization procedure to test the null hypothesis 
#     that there was no difference in the spatial distribution between years.
#     ...
#     we were testing only if the observed overlap was smaller than
#     random overlap, we considered this a one-tailed test. We
#     used a significance value of α=0.05 for one-way pairwise
#     comparisons. Following the null hypothesis above, p≥0.95
#     would indicate that the real observed value was at the other
#     end of the randomized distribution (on the right-end tail),
#     indicating significant clustering, or two species overlapping
#     more than expected by chance.
#     ... 
#     To do so, we calculated a p value as the proportion of randomized overlaps
#     that were smaller or equal to the mean observed overlap value 
#     for each pairwise comparison of species (null test)."
# Delord, K., et al. (2020). Movements of three alcid species breeding sympatrically
# in Saint Pierre and Miquelon, northwestern Atlantic Ocean. 
# Journal of Ornithology, 161(2), 359-371.
####

#### Each individual track is assigned winter year at random

## Get IDs of tracks and observed year (this will be shuffled, maintaing the 
## proportion of tracks per year)

df_years_sampled <- 
  tracks_sf %>% 
  st_drop_geometry() %>% 
  group_by(file_name) %>% 
  summarize(winter_year = unique(winter_year)) %>% 
  ungroup()

nrow(df_years_sampled)

## Create empty lists to store results
n_shuffles <- 1000
ll_overlap_years_UDOI_condF <- vector("list", length = n_shuffles)
ll_overlap_years_UDOI_95    <- vector("list", length = n_shuffles)
ll_overlap_years_UDOI_50    <- vector("list", length = n_shuffles)

Sys.time()

for (i in 1:1000) {
  
  ## Randomize the winter year
  df_years_sampled_new <-
    df_years_sampled %>% 
    mutate(winter_year_new = sample(winter_year, replace = FALSE))
  
  ## Assign to the dataset
  tracks_sf_laea_new <-
    tracks_sf %>% 
    left_join(df_years_sampled_new, by = c("file_name", "winter_year")) %>% 
    st_transform(crs = laea_proj)
  
  ## Convert to sp
  tracks_sp_laea_new <- as(tracks_sf_laea_new, "Spatial")
  
  ## Calculate kernels per year with randomized years
  KUD_YEAR <- kernelUD(tracks_sp_laea_new[, 11], h = 186000, same4all = TRUE, grid = 1000)
  
  ## UDOI Unconditional kernels overlap
  overlap_years_UDOI_condF <- kerneloverlaphr(KUD_YEAR, meth = "UDOI", conditional = FALSE)
  overlap_years_UDOI_condF[lower.tri(overlap_years_UDOI_condF)] <- NA
  diag(overlap_years_UDOI_condF) <- NA
  # ll_overlap_years_UDOI_condF[[i]] <- as.numeric(as.numeric(overlap_years_UDOI_condF) >= 0.9)
  # ll_overlap_years_UDOI_condF[[i]] <- as.numeric(overlap_years_UDOI_condF)
  ll_overlap_years_UDOI_condF[[i]] <- overlap_years_UDOI_condF
  
  ## UDOI 95% kernels overlap
  overlap_years_UDOI_95 <- kerneloverlaphr(KUD_YEAR, meth = "UDOI", conditional = TRUE, percent = 95)
  overlap_years_UDOI_95[lower.tri(overlap_years_UDOI_95)] <- NA
  diag(overlap_years_UDOI_95) <- NA
  # ll_ll_overlap_years_UDOI_95[[i]] <- as.numeric(as.numeric(overlap_years_UDOI_95) >= 0.9)
  # ll_overlap_years_UDOI_95[[i]] <- as.numeric(overlap_years_UDOI_95)
  ll_overlap_years_UDOI_95[[i]] <- overlap_years_UDOI_95
  
  ## UDOI 50% kernels overlap
  overlap_years_UDOI_50 <- kerneloverlaphr(KUD_YEAR, meth = "UDOI", conditional = TRUE, percent = 50)
  overlap_years_UDOI_50[lower.tri(overlap_years_UDOI_50)] <- NA
  diag(overlap_years_UDOI_50) <- NA
  # ll_ll_overlap_years_UDOI_50[[i]] <- as.numeric(as.numeric(overlap_years_UDOI_50) >= overlap_years_UDOI_50_obs)
  # ll_overlap_years_UDOI_50[[i]] <- as.numeric(overlap_years_UDOI_50)
  ll_overlap_years_UDOI_50[[i]] <- overlap_years_UDOI_50
  
  message(paste("Done ", i, "----"))
  message(Sys.time())
  
  if (i%%10 == 0) {
    Sys.sleep(5)
  }
  if (i%%100 == 0) {
    Sys.sleep(60*5)
    save(ll_overlap_years_UDOI_condF, file = paste0(my_path, "/RESULTS/ll_overlap_years_UDOI_condF.RData"))
    save(ll_overlap_years_UDOI_95, file = paste0(my_path, "/RESULTS/ll_overlap_years_UDOI_95.RData"))
    save(ll_overlap_years_UDOI_50, file = paste0(my_path, "/RESULTS/ll_overlap_years_UDOI_50.RData"))
    
  }
  
}


save(ll_overlap_years_UDOI_condF, file = paste0(my_path, "/RESULTS/ll_overlap_years_UDOI_condF.RData"))
save(ll_overlap_years_UDOI_95, file = paste0(my_path, "/RESULTS/ll_overlap_years_UDOI_95.RData"))
save(ll_overlap_years_UDOI_50, file = paste0(my_path, "/RESULTS/ll_overlap_years_UDOI_50.RData"))

Sys.time()


#### --------------------------------------------------------------------------~
#### CALCULATE THE OBSERVED VALUES OF OVERLAPS                              ----
#### --------------------------------------------------------------------------~

## Calculate kernel
KUD_YEAR_obs <- kernelUD(tracks_sp_laea[, 1], h = 186000, same4all = TRUE, grid = 1000)
image(KUD_YEAR_obs)

## UDOI Unconditional kernels overlap
overlap_years_UDOI_condF_obs <- kerneloverlaphr(KUD_YEAR_obs, meth = "UDOI", conditional = FALSE)
overlap_years_UDOI_condF_obs[lower.tri(overlap_years_UDOI_condF_obs)] <- NA
diag(overlap_years_UDOI_condF_obs) <- NA
# overlap_years_UDOI_condF_obs <- as.numeric(overlap_years_UDOI_condF_obs)
overlap_years_UDOI_condF_obs_df <- 
  data.frame("winter_years" = paste(rep(rownames(overlap_years_UDOI_condF_obs), times = 5), 
                                    rep(colnames(overlap_years_UDOI_condF_obs), each = 5), 
                                    sep = "_"), 
             "value_obs" = as.numeric(overlap_years_UDOI_condF_obs)) %>% 
  tidyr::drop_na()

## UDOI 95% kernels overlap
overlap_years_UDOI_95_obs <- kerneloverlaphr(KUD_YEAR_obs, meth = "UDOI", conditional = TRUE, percent = 95)
overlap_years_UDOI_95_obs[lower.tri(overlap_years_UDOI_95_obs)] <- NA
diag(overlap_years_UDOI_95_obs) <- NA
# overlap_years_UDOI_95_obs <- as.numeric(overlap_years_UDOI_95_obs)
overlap_years_UDOI_95_obs_df <- 
  data.frame("winter_years" = paste(rep(rownames(overlap_years_UDOI_95_obs), times = 5), 
                                    rep(colnames(overlap_years_UDOI_95_obs), each = 5), 
                                    sep = "_"), 
             "value_obs" = as.numeric(overlap_years_UDOI_95_obs)) %>% 
  tidyr::drop_na()


## UDOI 50% kernels overlap
overlap_years_UDOI_50_obs <- kerneloverlaphr(KUD_YEAR_obs, meth = "UDOI", conditional = TRUE, percent = 50)
overlap_years_UDOI_50_obs[lower.tri(overlap_years_UDOI_50_obs)] <- NA
diag(overlap_years_UDOI_50_obs) <- NA
# overlap_years_UDOI_50_obs <- as.numeric(overlap_years_UDOI_50_obs)
overlap_years_UDOI_50_obs_df <- 
  data.frame("winter_years" = paste(rep(rownames(overlap_years_UDOI_50_obs), times = 5), 
                                    rep(colnames(overlap_years_UDOI_50_obs), each = 5), 
                                    sep = "_"), 
             "value_obs" = as.numeric(overlap_years_UDOI_50_obs)) %>% 
  tidyr::drop_na()


#### --------------------------------------------------------------------------~
#### CALCULATE STATISTICS ----
#### --------------------------------------------------------------------------~

library(broom)
options(scipen = 999)
options(pillar.sigfig = 5)


#### --------------------------------------------------------------------------~
#### UDOI complete ----
#### --------------------------------------------------------------------------~

UDOI_condF_df <- 
  map_df(ll_overlap_years_UDOI_condF, as.data.frame) %>% 
  tibble::rownames_to_column(var = "winter_year_1") %>%
  dplyr::mutate(winter_year_1 = substr(winter_year_1, 1, 4)) %>% 
  tidyr::pivot_longer(!winter_year_1,
                      names_to = "winter_year_2",
                      values_to = "value") %>%
  tidyr::drop_na() %>%
  dplyr::mutate(winter_years = paste(winter_year_1, winter_year_2, sep = "_")) %>% 
  left_join(overlap_years_UDOI_condF_obs_df, by = "winter_years")


## t.test
UDOI_condF_ttest_less <-
  UDOI_condF_df %>% 
  split(as.character(.$winter_years)) %>%
  map(~ tidy(t.test(.x$value, mu = unique(.$value_obs), alternative = "less"))) %>% 
  bind_rows(.id = "winter_years")%>% 
  left_join(overlap_years_UDOI_condF_obs_df, by = "winter_years") %>% 
  dplyr::select(winter_years, value_obs, everything()) %>% 
  dplyr::mutate(p.value = format(round(p.value, 5), scientific=FALSE)) %>% 
  dplyr::select(winter_years, value_obs, estimate, p.value) %>% 
  mutate(overlap = "UDOI_all")

UDOI_condF_ttest_greater <-
  UDOI_condF_df %>%
  split(as.character(.$winter_years)) %>%
  map(~ tidy(t.test(.x$value, mu = unique(.$value_obs), alternative = "greater"))) %>%
  bind_rows(.id = "winter_years") %>%
  left_join(overlap_years_UDOI_condF_obs_df, by = "winter_years") %>%
  dplyr::select(winter_years, value_obs, everything()) %>%
  dplyr::mutate(p.value = format(round(p.value, 5), scientific=FALSE)) %>% 
  dplyr::select(winter_years, value_obs, estimate, p.value) %>% 
  mutate(overlap = "UDOI_all")



#### --------------------------------------------------------------------------~
#### KERNEL 95 ----
#### --------------------------------------------------------------------------~

UDOI_95_df <- 
  map_df(ll_overlap_years_UDOI_95, as.data.frame) %>% 
  tibble::rownames_to_column(var = "winter_year_1") %>%
  dplyr::mutate(winter_year_1 = substr(winter_year_1, 1, 4)) %>% 
  tidyr::pivot_longer(!winter_year_1,
                      names_to = "winter_year_2",
                      values_to = "value") %>%
  tidyr::drop_na() %>%
  dplyr::mutate(winter_years = paste(winter_year_1, winter_year_2, sep = "_")) %>% 
  left_join(overlap_years_UDOI_95_obs_df, by = "winter_years")

## t.test
UDOI_95_ttest_less <-
  UDOI_95_df %>% 
  split(as.character(.$winter_years)) %>%
  map(~ tidy(t.test(.x$value, mu = unique(.$value_obs), alternative = "less"))) %>% 
  bind_rows(.id = "winter_years")%>% 
  left_join(overlap_years_UDOI_95_obs_df, by = "winter_years") %>% 
  dplyr::select(winter_years, value_obs, everything()) %>% 
  dplyr::mutate(p.value = format(round(p.value, 5), scientific=FALSE)) %>% 
  dplyr::select(winter_years, value_obs, estimate, p.value) %>% 
  mutate(overlap = "UDOI_95")


UDOI_95_ttest_greater <-
  UDOI_95_df %>%
  split(as.character(.$winter_years)) %>%
  map(~ tidy(t.test(.x$value, mu = unique(.$value_obs), alternative = "greater"))) %>%
  bind_rows(.id = "winter_years") %>%
  left_join(overlap_years_UDOI_95_obs_df, by = "winter_years") %>%
  dplyr::select(winter_years, value_obs, everything()) %>%
  dplyr::mutate(p.value = format(round(p.value, 5), scientific=FALSE)) %>% 
  dplyr::select(winter_years, value_obs, estimate, p.value) %>% 
  mutate(overlap = "UDOI_95")



#### --------------------------------------------------------------------------~
#### KERNEL 50 ----
#### --------------------------------------------------------------------------~

UDOI_50_df <- 
  map_df(ll_overlap_years_UDOI_50, as.data.frame) %>% 
  tibble::rownames_to_column(var = "winter_year_1") %>%
  dplyr::mutate(winter_year_1 = substr(winter_year_1, 1, 4)) %>% 
  tidyr::pivot_longer(!winter_year_1,
                      names_to = "winter_year_2",
                      values_to = "value") %>%
  tidyr::drop_na() %>%
  dplyr::mutate(winter_years = paste(winter_year_1, winter_year_2, sep = "_")) %>% 
  left_join(overlap_years_UDOI_50_obs_df, by = "winter_years")

## t.test alternative hypothesis: true mean is less than observed value
## >> the observed is higher than expected from random
UDOI_50_ttest_less <-
  UDOI_50_df %>% 
  split(as.character(.$winter_years)) %>%
  map(~ tidy(t.test(.x$value, mu = unique(.$value_obs), alternative = "less"))) %>% 
  bind_rows(.id = "winter_years") %>% 
  left_join(overlap_years_UDOI_50_obs_df, by = "winter_years") %>% 
  # dplyr::select(winter_years, value_obs, everything() #%>% 
  dplyr::mutate(p.value = format(round(p.value, 5), scientific=FALSE)) %>% 
  dplyr::select(winter_years, value_obs, estimate, p.value) %>% 
  mutate(overlap = "UDOI_50")

UDOI_50_ttest_greater <-
  UDOI_50_df %>% 
  split(as.character(.$winter_years)) %>%
  map(~ tidy(t.test(.x$value, mu = unique(.$value_obs), alternative = "greater"))) %>% 
  bind_rows(.id = "winter_years")%>% 
  left_join(overlap_years_UDOI_50_obs_df, by = "winter_years") %>% 
  # dplyr::select(winter_years, value_obs, everything() #%>% 
  dplyr::mutate(p.value = format(round(p.value, 5), scientific=FALSE)) %>% 
  dplyr::select(winter_years, value_obs, estimate, p.value) %>% 
  mutate(overlap = "UDOI_50")


#### ONE TABLE OF ALL 
UDOI_ttests_less <-
  bind_rows(UDOI_50_ttest_less, 
            UDOI_95_ttest_less,
            UDOI_condF_ttest_less) %>% 
  dplyr::mutate(value_obs = round(value_obs, 3), 
                estimate = round(estimate, 3)) %>% 
  tidyr::pivot_wider(id_cols = "winter_years", 
                     names_from = overlap, 
                     values_from = c(value_obs, estimate, p.value)) %>% 
  dplyr::select(winter_years, ends_with("50"), ends_with("95"))


UDOI_ttests_greater <-
  bind_rows(UDOI_50_ttest_greater, 
            UDOI_95_ttest_greater,
            UDOI_condF_ttest_greater) %>% 
  dplyr::mutate(value_obs = round(value_obs, 3), 
                estimate = round(estimate, 3)) %>% 
  tidyr::pivot_wider(id_cols = "winter_years", 
                     names_from = overlap, 
                     values_from = c(value_obs, estimate, p.value)) %>% 
  dplyr::select(winter_years, ends_with("50"), ends_with("95"))



ft1 <- flextable::flextable(UDOI_ttests_less)
ft1 <- flextable::set_caption(ft1, "t-test mean less than observed")

ft2 <- flextable::flextable(UDOI_ttests_greater)
ft2 <- flextable::set_caption(ft2, "t-test mean greater than observed")

flextable::save_as_docx(ft1, ft2,
                        path = paste0(my_path, "/RESULTS/Table_overlaps_UDOI_final.docx"))


