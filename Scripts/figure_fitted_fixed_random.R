# Figure 7: Comparison of observed, predicted, fixed effects, random effects

library(tidyverse)
library(INLA)
library(sf)
library(sp)
library(cowplot)
library(rnaturalearth)


# Load R objects & data ----------

# Load model fit
fit_SR <- readRDS("Objects/SR_fit.rds")

# Load data stack used in model fitting
stack_SR <- readRDS("Objects/SR_stack.rds")

# Load data 
dat_SR <- read_csv("Data/SR_input.csv") %>%
  filter(Ntrials >= 25) %>%
  mutate(X_km = X_m/1000) %>% mutate(Y_km = Y_m/1000)

# Standardize covariates

# Function for standardization
stand_fun <- function(x){(x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)}

# Variables *not* to standardize
nostand <- c("damage_year", "X_m", "Y_m", "X_km", "Y_km", 
             "hotspot_bin", "hotspot_gt0", "Ntrials_gt0" )
nostand <- which(names(dat_SR) %in% nostand)

# Apply standardization to data frames
dat_SR <- cbind(dat_SR[,c(nostand)], apply(dat_SR[,-c(nostand)],2,stand_fun))


# Extract fitted values & linear predictors -------------

extract_fit_linpred <- function(fit, dat, stack) {

  # Observed data
  y = dat$hotspot_bin
  m = length(y)
  size = dat$Ntrials_gt0
  
  # Data indices
  data.index.z <- inla.stack.index(stack, "zobs")$data
  data.index.y <- inla.stack.index(stack, "yobs")$data
  
  # Vector of fitted values (success probabilities) for y>0 observations
  # Set to zero where y=0
  fitted_p = rep(0, m)
  fitted_p2 = fit$summary.fitted.values[data.index.y, "mean"][y > 0] 
  fitted_p[which(y > 0)] = fitted_p2
  
  # Fitted zero-probabilities
  fitted_p1 = fit$summary.fitted.values[data.index.z, "mean"]
  fitted_p0 = 1 - fitted_p1

  # Linear predictors
  linpred_p1 = fit$summary.linear.predictor[data.index.z, "mean"]
  linpred_p = rep(NA, m)
  linpred_p[which(y > 0)] = fit$summary.linear.predictor[data.index.y, "mean"][y > 0] 
  
  # Return data frame
  df <- data.frame("linpred_p1" = linpred_p1,
                   "linpred_p" = linpred_p,
                   "fitted_p1" = fitted_p1,
                   "fitted_p0" = fitted_p0,
                   "fitted_p" = fitted_p)
  return(df)
  
}

dat_SR <- cbind(dat_SR, extract_fit_linpred(fit_SR, dat_SR, stack_SR))


# Extract temporal random effects (random walks) ----------
z.rw1_year_SR <- fit_SR$summary.random$z.damage_year
u.rw1_year_SR <- fit_SR$summary.random$y.damage_year


# Extract spatio-temporal random fields -----------

# Extract coordinates
coords_SR <- select(dat_SR, X_km, Y_km) %>% distinct() %>% as.matrix()

# Build the boundary around the coordinates
boundary.loc_SR <- SpatialPoints(coords_SR)
bound_SR <- list(
  inla.nonconvex.hull(coordinates(boundary.loc_SR), 30),  # inner boundary
  inla.nonconvex.hull(coordinates(boundary.loc_SR), 70)) # outer boundary

# Build the SPDE mesh
mesh_SR <- inla.mesh.2d(boundary = bound_SR, 
                        cutoff = 5, max.edge = c(30, 200), min.angle=c(30, 21))

# Extract posterior mean of random fields
z.pm_SR <- fit_SR$summary.random$z$mean
u.pm_SR <- fit_SR$summary.random$u$mean

# Create projector matrix
projgrid_SR <- inla.mesh.projector(mesh_SR, coords_SR)


# Project random fields onto grid ---------

# Define data "groups" (in this case, years)
groups <- dat_SR$damage_year - (min(dat_SR$damage_year - 1))
ngroups <- length(unique(groups))

# Define the spatial fields
# This tells INLA which rows belong to the same year
field.z.idx <- inla.spde.make.index(name = 'z',
                                    n.spde = mesh_SR$n,
                                    n.group = ngroups)
field.u.idx <- inla.spde.make.index(name = 'u',
                                    n.spde = mesh_SR$n,
                                    n.group = ngroups)

# Loop through years (groups)
ncoords <- nrow(coords_SR)
years <- unique(dat_SR$damage_year)
nyears <- length(years)


z.pm.proj_SR <- data.frame("X_km" = rep(coords_SR[,1], times = nyears),
                           "Y_km" = rep(coords_SR[,2], times = nyears),
                           "damage_year" = rep(years, each = nrow(coords_SR)),
                           "z.pm" = NA)

u.pm.proj_SR <- data.frame("X_km" = rep(coords_SR[,1], times = nyears),
                           "Y_km" = rep(coords_SR[,2], times = nyears),
                           "damage_year" = rep(years, each = nrow(coords_SR)),
                           "u.pm" = NA)

for(i in 1:nyears){
  z.pm.proj_SR$z.pm[((i-1)*ncoords+1):(i*ncoords)] <- inla.mesh.project(projgrid_SR, z.pm_SR[field.z.idx$z.group==i])
  u.pm.proj_SR$u.pm[((i-1)*ncoords+1):(i*ncoords)] <- inla.mesh.project(projgrid_SR, u.pm_SR[field.u.idx$u.group==i])
}

# Join random effects to data frame
dat_SR <- dat_SR %>%
  left_join(z.pm.proj_SR) %>%
  left_join(u.pm.proj_SR) %>%
  left_join( mutate(z.rw1_year_SR, damage_year = ID) %>%
               mutate(z.rw = mean) %>%
               select(damage_year, z.rw) ) %>%
  left_join( mutate(u.rw1_year_SR, damage_year = ID) %>%
               mutate(u.rw = mean) %>%
               select(damage_year, u.rw) )


# Calculate fixed effects -------------------

# Calculate fixed effects by multiplying 
# coefficient estimates (post.mean) by covariates

# Coefficient estimates (post.mean, sorted)
all.coef <- fit_SR$summary.fixed %>%
  mutate(covariate = rownames(.)) %>%
  arrange(covariate)

# Covariates
all.covars <- dat_SR %>%
  mutate(tmin_norm_X_anom = tmin_norm * tmin_anom) %>%
  mutate(vpdmax_norm_X_anom = vpdmax_norm * vpdmax_anom) %>%
  select(AET_norm, HLI, Host_BA, Host_presence, Host_rich,
         tmin_norm, tmin_norm_X_anom, tmin_anom, 
         TWI, vpdmax_norm, vpdmax_norm_X_anom, vpdmax_anom)

# Multiply
fixed_p_multiply = as.numeric( as.matrix( select(all.covars, -Host_presence) ) %*% all.coef$mean[1:11] )
fixed_p1_multiply = as.numeric( as.matrix(all.covars) %*% all.coef$mean[12:23] )

dat_SR <- dat_SR %>%
  mutate(fixed_p_multiply = fixed_p_multiply) %>%
  mutate(fixed_p1_multiply = fixed_p1_multiply)


# Single year plotting of observed, predicted, and fixed versus random effects -------------

# CRS to use for plotting
my_crs <- CRS("ESRI:102039")

# Download US states vector data
usa_political <- ne_states(country = 'United States of America') %>% 
  st_as_sf() %>% st_transform(my_crs)


## OCCURRENCE, SOUTHERN ROCKIES ##

# get limits for plotting on logit scale
lim <- dat_SR %>% filter(damage_year == 2007) %>%
  mutate(randlim = z.rw + z.pm)
lim <- range(lim$randlim)

# get limits for SR region
xmin <- min(dat_SR$X_m)
xmax <- max(dat_SR$X_m)
ymin <- min(dat_SR$Y_m)
ymax <- max(dat_SR$Y_m)


p1 = ggplot() +
  geom_sf(data = usa_political, fill = NA, color="gray50", lwd = 0.1) +
  geom_raster(data = filter(dat_SR, damage_year == 2007),
              mapping = aes(X_m, Y_m, fill = factor(hotspot_bin))) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
  scale_fill_viridis_d(option = "rocket", name = "Observed\noccurrence ", end = 0.9) + 
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.height = unit(10, "pt"),
        legend.key.width = unit(10, "pt"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())

p2 = ggplot() +
  geom_sf(data = usa_political, fill = NA, color="gray50", lwd = 0.1) +
  geom_raster(filter(dat_SR, damage_year == 2007),
              mapping = aes(X_m, Y_m, fill = fitted_p1)) + 
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
  scale_fill_viridis_c(option = "rocket", name = "Predicted\np(occurrence) ", limits = c(0,1),
                       breaks = c(0,0.5,1), labels = c("0", "0.5", "1")) + 
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.height = unit(8, "pt"),
        legend.key.width = unit(8, "pt"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())

p3 = ggplot() +
  geom_sf(data = usa_political, fill = NA, color="gray50", lwd = 0.1) +
  geom_raster(filter(dat_SR, damage_year == 2007),
              mapping = aes(X_m, Y_m, fill = fixed_p1_multiply)) + 
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
  scale_fill_viridis_c(option = "rocket", limits = lim,
                       name = "Fixed effects \n(logit scale)") + 
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.height = unit(8, "pt"),
        legend.key.width = unit(8, "pt"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())

p4 = ggplot() +
  geom_sf(data = usa_political, fill = NA, color="gray50", lwd = 0.1) +
  geom_raster(filter(dat_SR, damage_year == 2007),
              mapping = aes(X_m, Y_m, fill = z.rw + z.pm)) + 
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax)) +
  scale_fill_viridis_c(option = "rocket", limits = lim,
                       name = "Random effects \n(logit scale)") + 
  theme_bw() +
  theme(legend.position = "bottom",
        legend.key.height = unit(8, "pt"),
        legend.key.width = unit(8, "pt"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())

plot_grid(p1, p2, p3, p4, nrow = 1, align = "h", label_x = 0.08, label_y = 0.95,
          labels = c("(a)", "(b)", "(c)", "(d)"), label_size = 10)

ggsave(paste0("Figures/obs_pred_fix_rand.png"), bg = "white",
       width = 6, height = 3.7, units = "in", dpi = 350)


