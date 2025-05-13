# Figure 1: study regions

# Code to install rnaturalearth
# devtools::install_github("ropenscilabs/rnaturalearthhires", force=TRUE)
# install.packages("rnaturalearthhires", repos="http://packages.ropensci.org")

library(rnaturalearth)
library(tidyverse)
library(sf)
library(raster)
require(ggspatial)

# CRS to use for plotting
my_crs <- CRS("ESRI:102039")

# Download US states vector data
usa_political <- ne_states(country = 'United States of America') %>% 
  st_as_sf() %>% st_transform(my_crs)

# Read in study region polygons
SR_region <- st_read("Data/Spatial/S_Rockies.shp") 
MR_region <- st_read("Data/Spatial/M_Rockies.shp")
C_region <- st_read("Data/Spatial/Cascades.shp")

# Read in hotspots data
SR_dat_full <- read_csv("Data/SR_input.csv") %>% drop_na(hotspot_bin)
MR_dat_full <- read_csv("Data/MR_input.csv") %>% drop_na(hotspot_bin)
C_dat_full <- read_csv("Data/C_input.csv") %>% drop_na(hotspot_bin)

SR_coords <- dplyr::select(SR_dat_full, X_m, Y_m) %>% distinct() %>%
  mutate(Region = "Southern Rockies")
MR_coords <- dplyr::select(MR_dat_full, X_m, Y_m) %>% distinct() %>%
  mutate(Region = "Middle Rockies")
C_coords <- dplyr::select(C_dat_full, X_m, Y_m) %>% distinct() %>%
  mutate(Region = "Cascades")

all_coords <- bind_rows(SR_coords, MR_coords, C_coords)

xmin <- min(all_coords$X_m)
xmax <- max(all_coords$X_m)
xdif <- xmax - xmin
ymin <- min(all_coords$Y_m)
ymax <- max(all_coords$Y_m)
ydif <- ymax - ymin


# Plot study area
ggplot() +
  geom_sf(data = usa_political, fill = NA, color= "gray50", lwd = 0.5) +
  geom_sf(data = SR_region, fill = NA, color = "#c36785", lwd = 0.5) +
  geom_sf(data = MR_region, fill = NA, color = "#7f64b9", lwd = 0.5) +
  geom_sf(data = C_region, fill = NA, color = "#a09344", lwd = 0.5) +
  geom_raster(data = SR_coords,
              aes(x = X_m, y = Y_m), fill = "#c36785", alpha = 0.5) + 
  geom_raster(data = MR_coords,
              aes(x = X_m, y = Y_m), fill = "#7f64b9", alpha = 0.5) + 
  geom_raster(data = C_coords,
              aes(x = X_m, y = Y_m), fill = "#a09344", alpha = 0.5) + 
  coord_sf(xlim = c(xmin - 150000, xmax + 100000), 
           ylim = c(ymin - 150000, ymax + 30000)) +
  annotate("text", x = xmin + 0.14*xdif, y = ymin + 0.72*ydif, 
           hjust = 0, color = "gray20", size = 3, fontface = "bold",
           label = "Cascades") +
  annotate("text", x = xmin + 0.75*xdif, y = ymin + 0.75*ydif, 
           hjust = 0, color = "gray20", size = 3, fontface = "bold",
           label = "Middle\nRockies") +
  annotate("text", x = xmin + 0.51*xdif, y = ymin + 0.21*ydif, 
           hjust = 0, color = "gray20", size = 3, fontface = "bold",
           label = "Southern\nRockies") +
  annotation_scale(location = "br", width_hint = 0.15) +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size = 8))


ggsave(paste0("Figures/study_regions.png"),
       width = 110, height = 110, units = "mm", dpi = 320)


