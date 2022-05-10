#Using tigris package to create a MN counties map
library(tigris)
library(ggplot2)
library(dplyr)


#allow tigris cache (speeds things up on later runs)
options(tigris_use_cache = TRUE)

#set sf (shapefile) option
options(tigris_class = "sf")

#grab 1:500k counties file (set with cb=T parameter)
#downloads from US Census Bureau
mn <- counties("MN", cb=T, progress_bar = FALSE)

# list of CIG counties
cig_cty <-
  c("Marshall",
    "Redwood",
    "Renville",
    "Clay",
    "Kittson",
    "Stearns",
    "Mower")

# flag the CIG counties
mn <-  mn %>% 
  mutate(INCLUDE = case_when(
    NAME %in% cig_cty ~ "Y", 
    TRUE ~ "N"
  ), 
  county_colors = case_when(
    NAME == "Marshall" ~ "RR",
    NAME == "Clay" ~ "RR",
    NAME == "Kittson" ~ "RR",
    NAME == "Renville" ~ "SW",
    NAME == "Redwood" ~ "SW",
    NAME == "Stearns" ~ "ST",
    NAME == "Mower" ~ "MW", 
    TRUE ~ "NOCIG"
  ))


#### MAP

# colors for county fill
# corresponds to scale color brewer qual palette 2
map_colors <- c(
  "RR" = "#377eb8",
  "SW" = "#984ea3",
  "ST" = "#e41a1c",
  "MW" = "#4daf4a", 
  "NOCIG" = "#f7f7f7"
)

# plotting
(cig_map <- ggplot(mn) +
  geom_sf(aes(fill = county_colors, geometry = geometry), show.legend = FALSE, color = "black") + 
  scale_fill_discrete(type = map_colors) +
  theme_void() + 
  theme(panel.grid.major = element_line(colour = 'white')))

# save 
 ggsave("mn_map_region_colors.png", width = 6, height = 6, units = "in")
