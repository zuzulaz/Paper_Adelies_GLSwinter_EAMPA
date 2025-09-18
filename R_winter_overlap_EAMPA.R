#### --------------------------------------------------------------------------~
#### WINTER DISTRIBUTION OF ADÉLIE PENGUINS                         ------------
#### OVERLAP OF LOCATIONS WITH THE PROPOSED EAMPA                   ------------
#### --------------------------------------------------------------------------~

#### Code for manuscript:
#### "Do penguins care about on-paper boundaries?
#### Conservation implications of spatio-temporal winter consistency 
#### in an Antarctic sentinel species: the Adélie penguin"

#### Zuzana Zajková
#### zajkova@icm.csic.es

#### --------------------------------------------------------------------------~
#### LIBRARIES ---- 
#### --------------------------------------------------------------------------~

library(tidyverse)
library(sf)

#### --------------------------------------------------------------------------~
#### 1) GETTING DATA ----
#### --------------------------------------------------------------------------~

#### --------------------------------------------------------------------------~
#### East Antarctic Marine Protected Area (EAMPA) proposal boundaries       ####
#### --------------------------------------------------------------------------~

mpa_ea_2015 <- st_read("/DATA_mpa/Earsmpa_2015_CCAMLR_XXXIV/Earsmpa_2015_CCAMLR_XXXIV.shp")

#### MacRobertson area
mpa_ea_2015_macr <- filter(mpa_ea_2015, mpa == 3)
#### Drygalski area
mpa_ea_2015_dryg <- filter(mpa_ea_2015, mpa == 5)
#### D’Urville Sea-Mertz area
mpa_ea_2015_mertz <- filter(mpa_ea_2015, mpa == 7)


#### --------------------------------------------------------------------------~
#### TRACKS                                                                 ####
#### --------------------------------------------------------------------------~

#### Data previously downloaded from Movebank and transformed to sf format
#### https://www.movebank.org/, study ID: 4778411162

load("Adelies_GLS_DDU_movebank.RData")

dim(tracks_sf)
# 23986     9


tracks_sf <- 
  tracks_sf %>% 
  dplyr::mutate(lon = as.numeric(longitude),
                lat = as.numeric(latitude), 
                month = lubridate::month(timestamp, label = FALSE), 
                month_label = lubridate::month(timestamp, label = TRUE, abbr = FALSE), 
                winter_year = as.character(lubridate::year(timestamp)), 
                gls_id = track_id,
                pit_id = bird_id,
                winter_year_gls_id = paste(winter_year, gls_id, sep ="_"), 
                file_name = paste("PYGADE", winter_year,bird_id, gls_id, sep ="_")) %>% 
  dplyr::select(winter_year, winter_year_gls_id, gls_id, 
                lon, lat, 
                date_gmt, time_gmt, month, month_label, 
                stage, file_name)

#### Reference data from Movebank
tracks_ref <- read_csv("Adelies_GLS_DDU_movebank_ref_table.csv",
                       col_types = cols(year_winter = col_character()))

#### Join
tracks_sf <- 
  tracks_sf %>% 
  left_join(tracks_ref, by = c("winter_year" = "year_winter", 
                               "gls_id" = "track_id"))



tracks_df <- 
  tracks_sf %>% 
  as.data.frame()


dim(tracks_sf)
dim(tracks_df)

length(unique(tracks_df$winter_year_gls_id))
# 61 (excluded winter_year_gls_id 3031)


tracks_sf_winter <- 
  tracks_sf %>% 
  dplyr::filter(stage == "wintering") 


tracks_sf_moult <- 
  tracks_sf %>% 
  dplyr::filter(stage == "moult") 




#### --------------------------------------------------------------------------~
#### Calculate the proportion of points within the MPA ----
#### --------------------------------------------------------------------------~

#### Calculated in longlat, Check the projections
st_crs(tracks_sf) == st_crs(mpa_ea_2015)


points_in_mertz <- 
  st_join(tracks_sf, mpa_ea_2015_mertz, join = st_within) %>% 
  dplyr::mutate(mpa = if_else(is.na(mpa),"mertz_out", "mertz_in"))

points_in_dryg <-
  st_join(tracks_sf, mpa_ea_2015_dryg, join = st_within) %>% 
  dplyr::mutate(mpa = if_else(is.na(mpa),"dryg_out", "dryg_in"))

points_in_macr <-
  st_join(tracks_sf, mpa_ea_2015_macr, join = st_within) %>% 
  dplyr::mutate(mpa = if_else(is.na(mpa),"macr_out", "macr_in"))



#### --------------------------------------------------------------------------~
#### Get proportions of points within MPA per winter_i and stage ----
#### --------------------------------------------------------------------------~

#### Mertz MPA
points_in_mertz_count <- 
  points_in_mertz %>% 
  as_tibble() %>% 
  group_by(winter_year, file_name, stage) %>% 
  count(mpa) %>% 
  tidyr::pivot_wider(names_from = c("mpa", "stage"), values_from = "n", values_fill = 0) %>% 
  dplyr::mutate(n_breeding_next = sum(mertz_in_breeding_next, mertz_out_breeding_next, na.rm = TRUE), 
                mertz_in_prop_breeding_next = mertz_in_breeding_next/n_breeding_next, 
                mertz_out_prop_breeding_next = mertz_out_breeding_next/n_breeding_next) %>% 
  dplyr::mutate(n_moult = sum(mertz_in_moult, mertz_out_moult, na.rm = TRUE), 
                mertz_in_prop_moult = mertz_in_moult/n_moult, 
                mertz_out_prop_moult = mertz_out_moult/n_moult) %>% 
  dplyr::mutate(n_premoult = sum(mertz_in_premoult, mertz_out_premoult, na.rm = TRUE), 
                mertz_in_prop_premoult = mertz_in_premoult/n_premoult, 
                mertz_out_prop_premoult = mertz_out_premoult/n_premoult) %>% 
  dplyr::mutate(n_wintering = sum(mertz_in_wintering, mertz_out_wintering, na.rm = TRUE), 
                mertz_in_prop_wintering = mertz_in_wintering/n_wintering, 
                mertz_out_prop_wintering = mertz_out_wintering/n_wintering) %>% 
  ungroup() %>% 
  mutate(pit_id = stringr::str_split(file_name, pattern = "_", simplify = TRUE)[,3])


#### Drygalski MPA
#### There are points in only during the wintering
points_in_dryg_count <- 
  points_in_dryg %>% 
  as_tibble() %>% 
  group_by(winter_year, file_name, stage) %>% 
  count(mpa) %>% 
  tidyr::pivot_wider(names_from = c("mpa", "stage"), values_from = "n", values_fill = 0) %>% 
  dplyr::mutate(n_wintering = sum(dryg_in_wintering, dryg_out_wintering, na.rm = TRUE), 
                dryg_in_prop_wintering = dryg_in_wintering/n_wintering, 
                dryg_out_prop_wintering = dryg_out_wintering/n_wintering) %>% 
  ungroup() %>% 
  mutate(pit_id = stringr::str_split(file_name, pattern = "_", simplify = TRUE)[,3])

#### Macr MPA
#### quick check
table(points_in_macr$mpa)
## All out.

#### Check map
ggplot() +
  geom_sf(data = points_in_mertz, aes(colour = mpa)) +
  geom_sf(data = filter(points_in_dryg, mpa == "dryg_in"), colour = "red") +
  geom_sf(data = mpa_ea_2015_mertz, colour = "black", fill = NA) +
  geom_sf(data = mpa_ea_2015_dryg, colour = "black", fill = NA)

#### --------------------------------------------------------------------------~
#### --------------------------------------------------------------------------~
#### Get stats  ----
#### --------------------------------------------------------------------------~
#### --------------------------------------------------------------------------~

library(lme4)


#### --------------------------------------------------------------------------~
#### MERTZ WINTER ---
#### --------------------------------------------------------------------------~
#### ONLY WHOLE TRACKS

points_in_mertz_count_58 <-
  points_in_mertz_count %>% 
  st_drop_geometry() %>% 
  dplyr::filter(!file_name %in% c("PYGADE_2015_972273000286176_3016",
                                  "PYGADE_2015_972273000285500_3023",
                                  "PYGADE_2015_972273000285376_3025")) 

length(table((points_in_mertz_count_58$file_name[points_in_mertz_count_58$mertz_in_prop_wintering > 0])))
length(table((points_in_mertz_count_58$file_name[points_in_mertz_count_58$mertz_in_prop_wintering == 0])))
points_in_mertz_count_58 %>% 
  filter(mertz_in_prop_wintering == 0) %>% 
  distinct(file_name)

points_in_mertz_count_58 %>% 
  pull(mertz_in_prop_wintering) %>% 
  summary() *100

points_in_mertz_count_58 %>% 
  pull(mertz_in_prop_wintering) %>% 
  sd() *100

points_in_mertz_count_58 %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ mean(.$mertz_in_prop_wintering*100))

points_in_mertz_count_58 %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ sd(.$mertz_in_prop_wintering*100))

points_in_mertz_count_58 %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ median(.$mertz_in_prop_wintering*100))


## only tracks with time > 0
points_in_mertz_count_58 %>% 
  filter(mertz_in_prop_wintering > 0) %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ mean(.$mertz_in_prop_wintering*100))

points_in_mertz_count_58 %>% 
  filter(mertz_in_prop_wintering > 0) %>% 
  pull(mertz_in_prop_wintering) %>% 
  summary() *100

## Plot
points_in_mertz_count_58 %>%
  # prepare data
  dplyr::select(winter_year, mertz_in_prop_wintering) %>%
  dplyr::group_by(winter_year) %>%
  dplyr::mutate(Mean = round(mean(mertz_in_prop_wintering*100), 1)) %>%
  dplyr::mutate(SD = round(sd(mertz_in_prop_wintering*100), 1)) %>%
  # start plot
  ggplot(aes(factor(winter_year), mertz_in_prop_wintering*100)) +
  geom_violin(trim=TRUE, color = "gray20")+ 
  geom_boxplot(width=0.1, fill="white", color = "gray20") +
  geom_text(aes(y=-2,label=paste("mean: ", Mean, sep = "")), size = 3, color = "black") +
  geom_text(aes(y=-5,label=paste("SD: ", SD, sep = "")), size = 3, color = "black") +
  theme(legend.position="none", 
        legend.title = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x = "Winter year", y = "mertz_in_prop_wintering*100")


#### Use Kruskal walis test, omiting repeated ID
kruskal.test(mertz_in_prop_wintering  ~ factor(winter_year), data = points_in_mertz_count_58)
# Kruskal-Wallis chi-squared = 8.1645, df = 4, p-value = 0.08573
# No sig.



#### --------------------------------------------------------------------------~
#### MERTZ MOULT ---
#### --------------------------------------------------------------------------~
#### ONLY WHOLE TRACKS, N=61

length(table((points_in_mertz_count$file_name[points_in_mertz_count$mertz_in_prop_moult > 0])))
length(table((points_in_mertz_count$file_name[points_in_mertz_count$mertz_in_prop_moult == 0])))
points_in_mertz_count %>% 
  filter(mertz_in_prop_moult == 0) %>% 
  distinct(file_name)

points_in_mertz_count %>% 
  filter(mertz_in_prop_moult > 0) %>% 
  distinct(file_name)

points_in_mertz_count %>% 
  pull(mertz_in_prop_moult) %>% 
  summary() *100

points_in_mertz_count %>% 
  pull(mertz_in_prop_moult) %>% 
  sd() *100

points_in_mertz_count %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ mean(.$mertz_in_prop_moult*100))

points_in_mertz_count %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ sd(.$mertz_in_prop_moult*100))

points_in_mertz_count %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ median(.$mertz_in_prop_moult*100))


## only tracks with time > 0
points_in_mertz_count %>% 
  filter(mertz_in_prop_moult > 0) %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ mean(.$mertz_in_prop_moult*100))

points_in_mertz_count %>% 
  filter(mertz_in_prop_moult > 0) %>% 
  pull(mertz_in_prop_moult) %>% 
  summary() *100

points_in_mertz_count %>% 
  filter(mertz_in_prop_moult > 0) %>% 
  split(.$file_name) %>% 
  purrr::map_dbl(~ min(.$n_moult))


## Plot
points_in_mertz_count %>%
  # prepare data
  dplyr::select(winter_year, mertz_in_prop_moult) %>%
  dplyr::group_by(winter_year) %>%
  dplyr::mutate(Mean = round(mean(mertz_in_prop_moult*100), 1)) %>%
  dplyr::mutate(SD = round(sd(mertz_in_prop_moult*100), 1)) %>%
  # start plot
  ggplot(aes(factor(winter_year), mertz_in_prop_moult*100)) +
  geom_violin(trim=TRUE, color = "gray20")+ 
  geom_boxplot(width=0.1, fill="white", color = "gray20") +
  geom_text(aes(y=-2,label=paste("mean: ", Mean, sep = "")), size = 3, color = "black") +
  geom_text(aes(y=-5,label=paste("SD: ", SD, sep = "")), size = 3, color = "black") +
  theme(legend.position="none", 
        legend.title = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x = "Winter year", y = "mertz_in_prop_moult*100")


#### Use Kruskal walis test, omiting repeated ID
kruskal.test(mertz_in_prop_moult ~ factor(winter_year), data = points_in_mertz_count)

#### Posthoc test
library(FSA)
dunn_result_moult <- dunnTest(mertz_in_prop_moult ~ factor(winter_year), 
                              data = points_in_mertz_count, 
                              method = "bonferroni")

# Extract the results
dunn_table_moult <- dunn_result_moult$res

# Export the Dunn's test results to Excel
library(writexl)
write_xlsx(dunn_table_moult, "/RESULTS/ADELIES_MPA_mertz_dunn_table_moult.xlsx")


#### --------------------------------------------------------------------------~
#### DRYGALSKI WINTER ---
#### --------------------------------------------------------------------------~

points_in_dryg_count_58 <-
  points_in_dryg_count %>% 
  st_drop_geometry() %>% 
  dplyr::filter(!file_name %in% c("PYGADE_2015_972273000286176_3016",
                                  "PYGADE_2015_972273000285500_3023",
                                  "PYGADE_2015_972273000285376_3025")) 

hist(points_in_dryg_count_58$dryg_in_prop_wintering*100, breaks = 100)
summary(points_in_dryg_count_58$dryg_in_prop_wintering)*100
sd(points_in_dryg_count_58$dryg_in_prop_wintering)*100

length(table((points_in_dryg_count_58$file_name[points_in_dryg_count_58$dryg_in_wintering > 0])))

points_in_dryg_count_58 %>% 
  st_drop_geometry() %>% 
  filter(dryg_in_wintering > 0) %>% 
  pull(dryg_in_prop_wintering) %>% 
  summary() *100

points_in_dryg_count_58 %>% 
  st_drop_geometry() %>% 
  filter(dryg_in_wintering > 0) %>% 
  pull(dryg_in_prop_wintering) %>% 
  sd() *100

points_in_dryg_count_58 %>% 
  st_drop_geometry() %>% 
  filter(dryg_in_wintering > 0) %>% 
  pull(dryg_in_prop_wintering) 

points_in_dryg_count_58 %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ mean(.$dryg_in_prop_wintering*100))

points_in_dryg_count_58 %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ sd(.$dryg_in_prop_wintering*100))

points_in_dryg_count_58 %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ median(.$dryg_in_prop_wintering*100))

## Plot
points_in_dryg_count_58 %>%
  # prepare data
  dplyr::select(winter_year, dryg_in_prop_wintering) %>%
  dplyr::group_by(winter_year) %>%
  dplyr::mutate(Mean = round(mean(dryg_in_prop_wintering*100), 1)) %>%
  dplyr::mutate(SD = round(sd(dryg_in_prop_wintering*100), 1)) %>%
  # start plot
  ggplot(aes(factor(winter_year), dryg_in_prop_wintering*100)) +
  geom_violin(trim=TRUE, color = "gray20")+ 
  geom_boxplot(width=0.1, fill="white", color = "gray20") +
  geom_text(aes(y=-2,label=paste("mean: ", Mean, sep = "")), size = 3, color = "black") +
  geom_text(aes(y=-5,label=paste("SD: ", SD, sep = "")), size = 3, color = "black") +
  theme(legend.position="none", 
        legend.title = element_blank(),
        panel.grid.minor = element_blank()) + 
  labs(x = "Winter year", y = "dryg_in_prop_wintering*100")


#### WINTER YEAR : Use Kruskal walis test, omiting repeated ID
kruskal.test(dryg_in_prop_wintering ~ factor(winter_year), data = points_in_dryg_count_58)
# NOT SIG.


#### --------------------------------------------------------------------------~
#### PER MONTHS ---
#### --------------------------------------------------------------------------~


#### Mertz MPA
points_in_mertz_count_month <- 
  points_in_mertz %>% 
  st_drop_geometry() %>% 
  as_tibble() %>% 
  filter(stage %in% c("moult", "wintering")) %>% 
  dplyr::filter(!(file_name %in% c("PYGADE_2015_972273000286176_3016",
                                   "PYGADE_2015_972273000285500_3023",
                                   "PYGADE_2015_972273000285376_3025") &
                    stage == "wintering")) %>%
  group_by(winter_year, file_name, month, month_label) %>%
  count(mpa) %>% 
  ungroup()

points_in_mertz_count_month_n <-
  points_in_mertz_count_month %>% 
  group_by(winter_year, file_name, month, month_label) %>% 
  summarize(n_month = sum(n)) %>% 
  ungroup() %>% 
  mutate(mpa_name = "mertz")

points_in_mertz_count_month_prop <- 
  points_in_mertz_count_month %>% 
  tidyr::pivot_wider(names_from = c("mpa"), values_from = "n", values_fill = 0) %>% 
  full_join(points_in_mertz_count_month_n, by = join_by(winter_year, file_name, month, month_label)) %>% 
  mutate(prop_in = mertz_in/n_month)

#### Winter months
points_in_mertz_count_month_winter <- 
  points_in_mertz %>% 
  st_drop_geometry() %>% 
  as_tibble() %>% 
  filter(stage %in% c("wintering")) %>% 
  dplyr::filter(!(file_name %in% c("PYGADE_2015_972273000286176_3016",
                                   "PYGADE_2015_972273000285500_3023",
                                   "PYGADE_2015_972273000285376_3025") &
                    stage == "wintering")) %>%
  # group_by(winter_year, file_name, month, month_label) %>% 
  group_by(file_name, month, month_label) %>%
  count(mpa) %>% 
  ungroup()

points_in_mertz_count_month_n_winter <-
  points_in_mertz_count_month_winter %>% 
  # group_by(winter_year, file_name, month, month_label) %>% 
  group_by(file_name, month, month_label) %>% 
  summarize(n_month = sum(n)) %>% 
  ungroup() %>% 
  mutate(mpa_name = "mertz")


points_in_mertz_count_month_prop_winter <- 
  points_in_mertz_count_month_winter %>% 
  tidyr::pivot_wider(names_from = c("mpa"), values_from = "n", values_fill = 0) %>% 
  full_join(points_in_mertz_count_month_n_winter, by = join_by(file_name, month, month_label)) %>% 
  mutate(prop_in = mertz_in/n_month) 




options(scipen=999)
n_distinct(points_in_mertz_count_month_prop_winter$file_name)


#### WINTER MONTHS
ggplot(data = points_in_mertz_count_month_prop_winter) +
  geom_boxplot(aes(x = factor(month), y = prop_in))

####Use Kruskal walis test, omiting repeated ID
kruskal.test(prop_in ~ factor(month), data = points_in_mertz_count_month_prop_winter)

#### Posthoc test
library(FSA)
dunn_result_winter_month <- dunnTest(prop_in ~ factor(month), 
                                     data = points_in_mertz_count_month_prop_winter, 
                                     method = "bonferroni")

# Extract the results
dunn_table_winter_month <- dunn_result_winter_month$res

# Export the Dunn's test results to Excel
library(writexl)
write_xlsx(dunn_table_winter_month, "/RESULTS/ADELIES_MPA_dunn_table_winter_months.xlsx")



#### Drygalski MPA
#### There are points in only during the wintering
points_in_dryg_count_month <- 
  points_in_dryg %>% 
  st_drop_geometry() %>% 
  as_tibble() %>% 
  filter(stage %in% c("moult", "wintering")) %>% 
  dplyr::filter(!(file_name %in% c("PYGADE_2015_972273000286176_3016",
                                   "PYGADE_2015_972273000285500_3023",
                                   "PYGADE_2015_972273000285376_3025") &
                    stage == "wintering")) %>%
  group_by(winter_year, file_name, month, month_label) %>% 
  count(mpa) 

points_in_dryg_count_month_n <-
  points_in_dryg_count_month %>% 
  group_by(winter_year, file_name, month, month_label) %>% 
  summarize(n_month = sum(n)) %>% 
  ungroup() %>% 
  mutate(mpa_name = "dryg")

points_in_dryg_count_month_prop <- 
  points_in_dryg_count_month %>% 
  tidyr::pivot_wider(names_from = c("mpa"), values_from = "n", values_fill = 0) %>% 
  full_join(points_in_dryg_count_month_n, by = join_by(winter_year, file_name, month, month_label)) %>% 
  mutate(prop_in = dryg_in/n_month)


points_in_mpa_count_month_prop <-
  bind_rows(points_in_dryg_count_month_prop, 
            points_in_mertz_count_month_prop) %>% 
  arrange(month) %>% 
  mutate(month_plot = forcats::as_factor(substr(month_label, 1,3)))


points_in_dryg_count_month_winter <- 
  points_in_dryg %>% 
  st_drop_geometry() %>% 
  as_tibble() %>% 
  filter(stage %in% c("wintering")) %>% 
  dplyr::filter(!(file_name %in% c("PYGADE_2015_972273000286176_3016",
                                   "PYGADE_2015_972273000285500_3023",
                                   "PYGADE_2015_972273000285376_3025") &
                    stage == "wintering")) %>%
  group_by(winter_year, file_name, month, month_label) %>% 
  count(mpa) 

points_in_dryg_count_month_n_winter <-
  points_in_dryg_count_month_winter %>% 
  group_by(winter_year, file_name, month, month_label) %>% 
  summarize(n_month = sum(n)) %>% 
  ungroup() %>% 
  mutate(mpa_name = "dryg")

points_in_dryg_count_month_prop_winter <- 
  points_in_dryg_count_month_winter %>% 
  tidyr::pivot_wider(names_from = c("mpa"), values_from = "n", values_fill = 0) %>% 
  full_join(points_in_dryg_count_month_n_winter, by = join_by(winter_year, file_name, month, month_label)) %>% 
  mutate(prop_in = dryg_in/n_month)


#### WINTER MONTHS
ggplot(data = points_in_dryg_count_month_prop_winter) +
  geom_boxplot(aes(x = factor(month), y = prop_in))

ggplot(data = points_in_dryg_count_month_prop_winter) +
  geom_boxplot(aes(x = factor(winter_year), y = prop_in))

ggplot(data = points_in_dryg_count_month_prop_winter) +
  geom_boxplot(aes(x = factor(month), y = prop_in, colour = winter_year)) +
  geom_point(aes(x = factor(month), y = prop_in, colour = winter_year)) +
  scale_colour_viridis_d(option = "magma", begin = 0.25, name = "Year") +
  scale_fill_viridis_d(option = "magma", begin = 0.25, name = "Year") 


points_in_dryg_count_month_prop_winter %>% 
  filter(prop_in > 0) %>% 
  split(.$month_label) %>% 
  map_dbl(~sd(.$prop_in))





#### Plot Figure S5 ----
points_in_mpa_count_month_prop <-
  bind_rows(points_in_dryg_count_month_prop, 
            points_in_mertz_count_month_prop) %>% 
  arrange(month) %>% 
  mutate(month_plot = forcats::as_factor(substr(month_label, 1,3)))

ggplot() +
  geom_boxplot(data = points_in_mpa_count_month_prop, 
               aes(x = factor(month), y = prop_in*100, 
                   fill = mpa_name, colour = mpa_name), 
               alpha = 0.6) +
  facet_wrap(~winter_year) +
  theme_bw()+
  scale_colour_manual(values = c("#E7C019", "#35899f"),  
                      labels = c("Drygalski", "D’Urville\nSea-Mertz"),
                      name = "Proposed MPA") +
  scale_fill_manual(values = c("#E7C019", "#3B9AB2"), 
                    labels = c("Drygalski", "D’Urville\nSea-Mertz"),
                    name = "Proposed MPA") +
  theme_bw()+
  theme(legend.position = c(0.88, 0.24),
        legend.direction = "vertical",
        legend.background = element_rect(fill = "#FFFFFF7E", colour = "#506A81", size = 0.15),
        legend.key.size = unit(20, "pt"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        axis.text = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "#506A81"), 
        strip.background = element_rect(fill = NA, colour = NA), 
        strip.text = element_text(size = 12)) +
  labs(x = "Month", y = "Proportion of locations within (%)")

ggsave(plot = last_plot(), 
       filename = paste0(my_path, "/PLOTS/MAPS/Adelies_Fig_4_prop_mpa_lowres.png"), 
       width = 10, height = 6, units = "in", dpi = 300)

#### --------------------------------------------------------------------------~
#### Table 3 Proportion of locations within EAMPA per stage ----
#### --------------------------------------------------------------------------~

#### Table Mertz ----
table_mertz_all_moult <-
  points_in_mertz %>% 
  dplyr::filter(stage %in% c("moult")) %>% 
  as_tibble() %>% 
  group_by(file_name, stage) %>% 
  count(mpa) %>% 
  tidyr::pivot_wider(names_from = c("mpa"), values_from = "n", values_fill = 0) %>% 
  ungroup() %>% 
  dplyr::mutate(total = mertz_in + mertz_out, 
                mertz_in_prop = mertz_in / total, 
                pmpa = "mertz") %>% 
  group_by(stage) %>% 
  dplyr::summarise(mertz_in_prop_mean = round(mean(mertz_in_prop*100, na.rm = TRUE), 1),
                   mertz_in_prop_sd = round(sd(mertz_in_prop*100, na.rm = TRUE), 1),
                   mertz_in_prop_n = length(mertz_in_prop)) %>% 
  ungroup() %>% 
  dplyr::mutate(mertz_in_prop_msd = paste0(mertz_in_prop_mean, " ± ", 
                                           mertz_in_prop_sd)) %>% 
  dplyr::select(-c(mertz_in_prop_mean, mertz_in_prop_sd)) %>% 
  dplyr::mutate(winter_year = "All") %>% 
  tidyr::pivot_wider(names_from = c("stage"), 
                     values_from = c("mertz_in_prop_msd", "mertz_in_prop_n"), 
                     values_fill = NA) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column() 

table_mertz_all_winter <-
  points_in_mertz %>% 
  dplyr::filter(stage %in% c("wintering")) %>% 
  dplyr::filter(!(file_name %in% c("PYGADE_2015_972273000286176_3016",
                                   "PYGADE_2015_972273000285500_3023",
                                   "PYGADE_2015_972273000285376_3025"))) %>%
  as_tibble() %>% 
  group_by(file_name, stage) %>% 
  count(mpa) %>% 
  tidyr::pivot_wider(names_from = c("mpa"), values_from = "n", values_fill = 0) %>% 
  ungroup() %>% 
  dplyr::mutate(total = mertz_in + mertz_out, 
                mertz_in_prop = mertz_in / total, 
                pmpa = "mertz") %>% 
  group_by(stage) %>% 
  dplyr::summarise(mertz_in_prop_mean = round(mean(mertz_in_prop*100, na.rm = TRUE), 1),
                   mertz_in_prop_sd = round(sd(mertz_in_prop*100, na.rm = TRUE), 1),
                   mertz_in_prop_n = length(mertz_in_prop)) %>% 
  ungroup() %>% 
  dplyr::mutate(mertz_in_prop_msd = paste0(mertz_in_prop_mean, " ± ", 
                                           mertz_in_prop_sd)) %>% 
  dplyr::select(-c(mertz_in_prop_mean, mertz_in_prop_sd)) %>% 
  dplyr::mutate(winter_year = "All")%>% 
  tidyr::pivot_wider(names_from = c("stage"), 
                     values_from = c("mertz_in_prop_msd", "mertz_in_prop_n"), 
                     values_fill = NA) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column() 
#### Change column names and remove the row
colnames(table_mertz_all_winter) <- as.character(table_mertz_all_winter[1, ])
table_mertz_all_winter <- table_mertz_all_winter[-1, ]
colnames(table_mertz_all_moult) <- as.character(table_mertz_all_moult[1, ])
table_mertz_all_moult <- table_mertz_all_moult[-1, ]


table_mertz_all <- 
  bind_rows(table_mertz_all_moult, 
            table_mertz_all_winter)

table_mertz_moult <-
  points_in_mertz %>% 
  as_tibble() %>% 
  group_by(winter_year, file_name, stage) %>% 
  count(mpa) %>% 
  tidyr::pivot_wider(names_from = c("mpa"), values_from = "n", values_fill = 0) %>% 
  ungroup() %>% 
  dplyr::filter(stage %in% c("moult")) %>% 
  dplyr::mutate(total = mertz_in + mertz_out, 
                mertz_in_prop = mertz_in / total, 
                pmpa = "mertz") %>% 
  group_by(winter_year, stage) %>% 
  dplyr::summarise(mertz_in_prop_mean = round(mean(mertz_in_prop*100, na.rm = TRUE), 1),
                   mertz_in_prop_sd = round(sd(mertz_in_prop*100, na.rm = TRUE), 1),
                   mertz_in_prop_n = length(mertz_in_prop)) %>% 
  ungroup() %>% 
  dplyr::mutate(mertz_in_prop_msd = paste0(mertz_in_prop_mean, " ± ", 
                                           mertz_in_prop_sd)) %>% 
  dplyr::select(-c(mertz_in_prop_mean, mertz_in_prop_sd)) %>% 
  # t()
  tidyr::pivot_wider(names_from = c("stage"), 
                     values_from = c("mertz_in_prop_msd", "mertz_in_prop_n"), 
                     values_fill = NA) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column()
#### Change column names and remove the row
colnames(table_mertz_moult) <- as.character(table_mertz_moult[1, ])
table_mertz_moult <- table_mertz_moult[-1, ]


table_mertz_winter <-
  points_in_mertz %>% 
  dplyr::filter(!(file_name %in% c("PYGADE_2015_972273000286176_3016",
                                   "PYGADE_2015_972273000285500_3023",
                                   "PYGADE_2015_972273000285376_3025"))) %>%
  as_tibble() %>% 
  group_by(winter_year, file_name, stage) %>% 
  count(mpa) %>% 
  tidyr::pivot_wider(names_from = c("mpa"), values_from = "n", values_fill = 0) %>% 
  ungroup() %>% 
  dplyr::filter(stage %in% c("wintering")) %>% 
  dplyr::mutate(total = mertz_in + mertz_out, 
                mertz_in_prop = mertz_in / total, 
                pmpa = "mertz") %>% 
  group_by(winter_year, stage) %>% 
  dplyr::summarise(mertz_in_prop_mean = round(mean(mertz_in_prop*100, na.rm = TRUE), 1),
                   mertz_in_prop_sd = round(sd(mertz_in_prop*100, na.rm = TRUE), 1),
                   mertz_in_prop_n = length(mertz_in_prop)) %>% 
  ungroup() %>% 
  dplyr::mutate(mertz_in_prop_msd = paste0(mertz_in_prop_mean, " ± ", 
                                           mertz_in_prop_sd)) %>% 
  dplyr::select(-c(mertz_in_prop_mean, mertz_in_prop_sd)) %>% 
  # t()
  tidyr::pivot_wider(names_from = c("stage"), 
                     values_from = c("mertz_in_prop_msd", "mertz_in_prop_n"), 
                     values_fill = NA) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column()

#### Change column names and remove the row
colnames(table_mertz_winter) <- as.character(table_mertz_winter[1, ])
table_mertz_winter <- table_mertz_winter[-1, ]


table_mertz_all <- 
  bind_rows(table_mertz_all_moult, 
            table_mertz_all_winter)

table_mertz <-
  bind_rows(table_mertz_moult, 
            table_mertz_winter) %>% 
  left_join(table_mertz_all, by = join_by(winter_year))



#### Table Drygalski ----
table_dryg_all_winter <-
  points_in_dryg %>% 
  dplyr::filter(stage %in% c("wintering")) %>% 
  dplyr::filter(!(file_name %in% c("PYGADE_2015_972273000286176_3016",
                                   "PYGADE_2015_972273000285500_3023",
                                   "PYGADE_2015_972273000285376_3025"))) %>%
  as_tibble() %>% 
  group_by(file_name, stage) %>% 
  count(mpa) %>% 
  tidyr::pivot_wider(names_from = c("mpa"), values_from = "n", values_fill = 0) %>% 
  ungroup() %>% 
  dplyr::mutate(total = dryg_in + dryg_out, 
                dryg_in_prop = dryg_in / total, 
                pmpa = "dryg") %>% 
  group_by(stage) %>% 
  dplyr::summarise(dryg_in_prop_mean = round(mean(dryg_in_prop*100, na.rm = TRUE), 1),
                   dryg_in_prop_sd = round(sd(dryg_in_prop*100, na.rm = TRUE), 1),
                   dryg_in_prop_n = length(dryg_in_prop)) %>% 
  ungroup() %>% 
  dplyr::mutate(dryg_in_prop_msd = paste0(dryg_in_prop_mean, " ± ", 
                                          dryg_in_prop_sd)) %>% 
  dplyr::select(-c(dryg_in_prop_mean, dryg_in_prop_sd)) %>% 
  dplyr::mutate(winter_year = "All")%>% 
  tidyr::pivot_wider(names_from = c("stage"), 
                     values_from = c("dryg_in_prop_msd", "dryg_in_prop_n"), 
                     values_fill = NA) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column() 
#### Change column names and remove the row
colnames(table_dryg_all_winter) <- as.character(table_dryg_all_winter[1, ])
table_dryg_all_winter <- table_dryg_all_winter[-1, ]


table_dryg_winter <-
  points_in_dryg %>% 
  dplyr::filter(!(file_name %in% c("PYGADE_2015_972273000286176_3016",
                                   "PYGADE_2015_972273000285500_3023",
                                   "PYGADE_2015_972273000285376_3025"))) %>%
  as_tibble() %>% 
  group_by(winter_year, file_name, stage) %>% 
  count(mpa) %>% 
  tidyr::pivot_wider(names_from = c("mpa"), values_from = "n", values_fill = 0) %>% 
  ungroup() %>% 
  dplyr::filter(stage %in% c("wintering")) %>% 
  dplyr::mutate(total = dryg_in + dryg_out, 
                dryg_in_prop = dryg_in / total, 
                pmpa = "dryg") %>% 
  group_by(winter_year, stage) %>% 
  dplyr::summarise(dryg_in_prop_mean = round(mean(dryg_in_prop*100, na.rm = TRUE), 1),
                   dryg_in_prop_sd = round(sd(dryg_in_prop*100, na.rm = TRUE), 1),
                   dryg_in_prop_n = length(dryg_in_prop)) %>% 
  ungroup() %>% 
  dplyr::mutate(dryg_in_prop_msd = paste0(dryg_in_prop_mean, " ± ", 
                                          dryg_in_prop_sd)) %>% 
  dplyr::select(-c(dryg_in_prop_mean, dryg_in_prop_sd)) %>% 
  # t()
  tidyr::pivot_wider(names_from = c("stage"), 
                     values_from = c("dryg_in_prop_msd", "dryg_in_prop_n"), 
                     values_fill = NA) %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column()

#### Change column names and remove the row
colnames(table_dryg_winter) <- as.character(table_dryg_winter[1, ])
table_dryg_winter <- table_dryg_winter[-1, ]


table_dryg <-
  table_dryg_winter %>% 
  left_join(table_dryg_all_winter, by = join_by(winter_year))


#### Join both ----
table_prop_mpa_3 <-
  bind_rows(table_mertz, 
            table_dryg)

#### Save
save(table_prop_mpa_3, 
     file = "/RESULTS/ADELIES_table_3_prop_mpa.RData")

write.csv(table_prop_mpa_3, 
          file = "/RESULTS/ADELIES_table_3_prop_mpa.csv", 
          row.names = FALSE)







#### INDIVIDUAL TABLE FOR SUPPL ----

table_mertz_ind <-
  points_in_mertz_count %>% 
  st_drop_geometry() %>% 
  dplyr::select(winter_year, file_name, pit_id,
                mertz_in_prop_moult,
                mertz_in_prop_wintering) %>% 
  dplyr::mutate(mertz_in_prop_moult = round(mertz_in_prop_moult*100, 1),
                mertz_in_prop_wintering = round(mertz_in_prop_wintering*100, 1))


table_dryg_ind <-
  points_in_dryg_count %>% 
  st_drop_geometry() %>% 
  dplyr::select(winter_year, file_name, pit_id,
                dryg_in_prop_wintering)


table_prop_ind <-
  table_mertz_ind %>%
  left_join(table_dryg_ind, by = join_by(winter_year, file_name, pit_id)) %>% 
  group_by(file_name) %>%
  dplyr::mutate(gls_id = unlist(strsplit(file_name, "_"))[4],
                winter_year_gls_id = paste(winter_year, gls_id, sep ="_")) %>% 
  ungroup() %>% 
  ## Create pseudo ID (to avoid long pit_id) %>% 
  dplyr::arrange(winter_year_gls_id) %>% 
  dplyr::mutate(pseudo_id = paste0("ID_", as.numeric(forcats::fct_inorder(pit_id))))


#### Save
save(table_prop_ind, 
     file = "/RESULTS/ADELIES_table_prop_mpa_stage_ind.RData")

write.csv(table_prop_ind, 
          file = "/RESULTS/ADELIES_table_prop_mpa_stage_ind.csv", 
          row.names = FALSE)


#### TABLE OF GENERAL USE OF EAMPA ----
#### ONLY FOR WINTER

count_mertz_ind <-
  points_in_mertz_count %>% 
  st_drop_geometry() %>% 
  dplyr::select(winter_year, file_name, pit_id,
                mertz_in_moult, mertz_out_moult, n_moult,
                mertz_in_wintering, mertz_out_wintering, n_wintering)


count_dryg_ind <-
  points_in_dryg_count %>% 
  st_drop_geometry() %>% 
  dplyr::select(winter_year, file_name, pit_id,
                dryg_in_wintering, dryg_out_wintering, n_wintering)


count_prop_ind_eampa <-
  count_mertz_ind %>%
  left_join(count_dryg_ind, by = join_by(winter_year, file_name, pit_id, n_wintering)) %>% 
  group_by(file_name) %>%
  dplyr::mutate(gls_id = unlist(strsplit(file_name, "_"))[4],
                winter_year_gls_id = paste(winter_year, gls_id, sep ="_")) %>% 
  ungroup() %>% 
  #### Calculations 
  dplyr::mutate(eampa_in_winter = (mertz_in_wintering + dryg_in_wintering), 
                eampa_in_prop_winter = eampa_in_winter/n_wintering)

count_prop_ind_eampa_58 <-
  count_prop_ind_eampa %>% 
  dplyr::filter(!(file_name %in% c("PYGADE_2015_972273000286176_3016",
                                   "PYGADE_2015_972273000285500_3023",
                                   "PYGADE_2015_972273000285376_3025"))) 


count_prop_ind_eampa_58 %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ median(.$eampa_in_prop_winter*100))


count_prop_ind_eampa_58 %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ mean(.$eampa_in_prop_winter*100))


count_prop_ind_eampa_58 %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ sd(.$eampa_in_prop_winter*100))

count_prop_ind_eampa_58 %>% 
  split(.$winter_year) %>% 
  purrr::map_dbl(~ length(.$eampa_in_prop_winter*100))

summary(count_prop_ind_eampa_58$eampa_in_prop_winter*100)
sd(count_prop_ind_eampa_58$eampa_in_prop_winter*100)
range(count_prop_ind_eampa_58$eampa_in_prop_winter*100)
hist(count_prop_ind_eampa_58$eampa_in_prop_winter*100, breaks = 20)
