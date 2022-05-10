
# setup -------------------------------------------------------------------

library(tidyverse)
library(glue)


# Climate data ------------------------------------------------------------
### citation: NOAA National Centers for Environmental information, 
# Climate at a Glance: Divisional Time Series, published April 2022, 
# retrieved on May 3, 2022 from https://www.ncdc.noaa.gov/cag/

#### Precipitation 
sw.map <- read_csv("MAP_mn_div7_southwest_redwood_renville_1991_2020.csv", skip = 4) %>% 
  mutate(division = "MN_div_7", 
         region = "SW")

mw.map <- read_csv("MAP_mn_div9_southeast_mower_1991_2020.csv", skip = 4) %>% 
  mutate(division = "MN_div_9", 
         region = "SE")

st.map <- read_csv("MAP_mn_div5_central_stearns_1991_2020.csv", skip = 4) %>% 
  mutate(division = "MN_div_5", 
         region = "ST")

rr.map <- read_csv("MAP_mn_div1_northwest_rrv_1991_2020.csv", skip = 4) %>% 
  mutate(division = "MN_div_1", 
         region = "RR")

map.all <- bind_rows(sw.map, mw.map, st.map, rr.map) %>% 
  rename(precipitation_in = Value)


#### Temperature 
sw.mat <- read_csv("MAT_mn_div7_southwest_redwood_renville_1991_2020.csv", skip = 4) %>% 
  mutate(division = "MN_div_7", 
         region = "SW")

mw.mat <- read_csv("MAT_mn_div9_southeast_mower_1991_2020.csv", skip = 4) %>% 
  mutate(division = "MN_div_9", 
         region = "SE")

st.mat <- read_csv("MAT_mn_div5_central_stearns_1991_2020.csv", skip = 4) %>% 
  mutate(division = "MN_div_5", 
         region = "ST")

rr.mat <- read_csv("MAT_mn_div1_northwest_rrv_1991_2020.csv", skip = 4) %>% 
  mutate(division = "MN_div_1", 
         region = "RR")

mat.all <- bind_rows(sw.mat, mw.mat, st.mat, rr.mat) %>% 
  rename(temperature_fahr = Value)



# Convert to metric -------------------------------------------------------

# precip inches to mm
map.all$precip.mm <- map.all$precipitation_in * 25.4

# temp fahrenheit to celsius
mat.all$temperature_celsius <- (mat.all$temperature_fahr - 32) * (5/9)


# calculate 30 year normals for MAP and MAT -------------------------------
# 1991-2020
map.norms <- map.all %>% 
  group_by(region) %>% 
  summarise(mean.annual.precip.mm = mean(precip.mm, na.rm = T))

mat.norms <- mat.all %>% 
  group_by(region) %>% 
  summarise(mean.annual.temp.c = mean(temperature_celsius))

# join together
norms.all <- inner_join(map.norms, mat.norms, by = "region")

# save 
write_csv(norms.all, "map_mat_normals_cig_regions.csv") 


# SOM, pH, and clay data --------------------------------------------------

lab.dat <- read_csv("../cig-main/cig_lab_data_all_20220408.csv")

# keeping 2019 samples here b/c that aligns w/ PLFA and EEA data presented
spc <- lab.dat %>% 
  select(sample_id, region, site, treatment, position, clay4hr_perc, pH, OM) %>% 
  rename(clay_perc = clay4hr_perc,
         som_loi = OM) %>% 
  filter(treatment == "SH" |
           treatment == "CV", 
         sample_id %in% as.character(c(1:243)))

# calculate means and sds 
spc.summary <- spc %>%
  group_by(region) %>%
  summarise(across(
    .cols = c(clay_perc, pH, som_loi),
    .fns = list(mean = mean, sd = sd), 
    .names = "{.col}.{.fn}"
  ))

# save 

write_csv(spc.summary, "regional_mean_sd_ph_som_clay.csv")


