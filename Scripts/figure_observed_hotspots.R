# Figures 2, 3, 4: observed hotspots

# Code to install rnaturalearth
# devtools::install_github("ropenscilabs/rnaturalearthhires", force=TRUE)
# install.packages("rnaturalearthhires", repos="http://packages.ropensci.org")

library(rnaturalearth)
library(tidyverse)
library(sf)
library(sp)
library(raster)
library(fasterize)

# Specify region
# 1 = Southern Rockies
# 2 = Cascades
# 3 = Middle Rockies
region = 1
prefix <- c("SR", "C", "MR")[region]

# CRS to use for plotting
my_crs <- CRS("ESRI:102039")

# Download US states vector data
usa_political <- ne_states(country = 'United States of America') %>% 
  st_as_sf() %>% st_transform(my_crs)

# Read in data
dat_full <- read_csv(paste0("Data/", prefix, "_input.csv")) %>%
  drop_na(hotspot_bin)
years <- sort(unique(dat_full$damage_year))
coords <- dplyr::select(dat_full, X_m, Y_m) %>% distinct()

xmin <- min(dat_full$X_m)
xmax <- max(dat_full$X_m)
xdif <- xmax - xmin
ymin <- min(dat_full$Y_m)
ymax <- max(dat_full$Y_m)
ydif <- ymax - ymin


# Plot observed spatio-temporal patterns ----

ggplot() +
  geom_sf(data = usa_political, fill = NA, color="gray50", lwd = 0.1) +
  geom_raster(data = filter(dat_full, hotspot_bin == 0),
              aes(x = X_m, y = Y_m), fill = "gray80") +
  geom_raster(data = filter(dat_full, hotspot_bin == 1),
              aes(x = X_m, y = Y_m,
                  fill = log10(hotspot_gt0 / Ntrials_gt0))) + 
  scale_fill_viridis_c(name = "Hotspot prevalence",
                       option = "magma",
                       limits = c(-2,0),
                       labels = c("0.01", "0.03", "0.1", "0.3", "1")) + 
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
  facet_wrap(~damage_year, nrow = ifelse(prefix == "MR", 4, 3)) +
  theme_bw() +
  theme(legend.position = c(0.7, -0.1),
        legend.title = element_text(size = 8, vjust = 1),
        legend.key.height = unit(8, "pt"),
        legend.key.width = unit(12, "pt"),
        legend.direction = "horizontal",
        legend.text = element_text(size = 7, angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(size = 7, margin = margin(t = 0, b = 1)),
        strip.background = element_blank(),
        panel.spacing = unit(2, "pt"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = margin(t = 5.5, r = 5.5, l = 5.5,
                             b = 50))

ggsave(paste0("Figures/hotspots_spatiotemporal_", prefix, ".png"),
       width = 4.25, height = 5.5, units = "in", dpi = 350)


# Plot observed temporal patterns ----

# plot widths vary by region
wid = ifelse(prefix == "SR", 3.1, ifelse(prefix == "MR", 3.5, 3.3))

# Hotspot prevalence

ggplot(dat_full, aes(x = as.factor(damage_year), y = hotspot_gt0 / Ntrials_gt0,
                     color = log10(hotspot_gt0 / Ntrials_gt0))) + 
  geom_jitter(alpha=0.05, width=0.15) +
  geom_boxplot(width=0.7, outlier.shape = NA, fill = NA) +
  labs(x = "Year", y = "Prevalence\n(proportion of\n510-m subcells)") +
  scale_y_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1),
                labels = c("0.01", "0.03", "0.1", "0.3", "1")) +
  scale_color_viridis_c(name = "Prevalence\n(proportion of\n510-m subcells)",
                        option = "magma",
                        limits = c(-2,0),
                        labels = c("0.01", "0.03", "0.1", "0.3", "1")) + 
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7, angle = 135, vjust = 0, hjust = 0),
        axis.title.y = element_text(size = 8),
        plot.margin = margin(t = 8, r = 5.5, l = 5.5, b = 5.5))

ggsave(paste0("Figures/hotspots_temporal_", prefix, ".png"),
       width = wid, height = 1, units = "in", dpi = 350)


# Hotspot occurrence

dat_stacked_bin <- dat_full %>%
  group_by(damage_year) %>%
  summarize(total_n = length(hotspot_bin),
            hotspot_bin_1 = sum(hotspot_bin)) %>%
  mutate(hotspot_bin_0 = total_n - hotspot_bin_1) %>%
  pivot_longer(cols = starts_with("hotspot_bin_"),
               names_to = "hotspot_bin",
               names_prefix = "hotspot_bin_",
               values_to = "count") %>%
  mutate(prop = count / total_n) %>%
  mutate(area_km2 = count * 5.1 * 5.1)

ggplot(filter(dat_stacked_bin, hotspot_bin == 1), 
       aes(x = as.factor(damage_year), y = prop), color = "gray20") + 
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0,1)) +
  labs(y = "Occurrence\n(proportion of\n 5.1-km cells)") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 7, angle = 135, vjust = 0, hjust = 0),
        axis.title.y = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 8, r = 5.5, l = 5.5, b = 5.5))

ggsave(paste0("Figures/hotspots_temporal_binary_", prefix, ".png"),
       width = wid, height = 1.25, units = "in", dpi = 350)



