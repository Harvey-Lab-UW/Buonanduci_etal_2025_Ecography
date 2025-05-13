## Spatio-temporal model-fitting with INLA

library(tidyverse)
library(INLA)
library(sp)

## Read in data -------

# Specify region
# 1 = Southern Rockies
# 2 = Cascades
# 3 = Middle Rockies
region = 1
prefix <- c("SR", "C", "MR")[region]

dat <- read_csv(paste0("Data/", prefix, "_input.csv"))

# Filter to only include cells at least 25% (a) surveyed and (b) with host cover
dat <- dat %>%
  filter(Ntrials >= 25)

# Rescale coordinates to kilometers
dat <- mutate(dat, X_km = X_m/1000) %>%
  mutate(Y_km = Y_m/1000)

# Standardize covariates

# Function for standardization
stand_fun <- function(x){(x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)}

# Variables *not* to standardize
nostand <- c("damage_year", "X_m", "Y_m", "X_km", "Y_km", 
             "hotspot_bin", "hotspot_gt0", "Ntrials_gt0" )
nostand <- which(names(dat) %in% nostand)

# Apply standardization to data frame
dat_stand <- cbind(dat[,c(nostand)], apply(dat[,-c(nostand)],2,stand_fun))


## Make SPDE mesh and SPDE object  -------

# Extract coordinates
coords <- select(dat_stand, X_km, Y_km) %>% distinct() %>% as.matrix()

# Build the boundary around the coordinates
boundary.loc <- SpatialPoints(coords)
bound <- list(
  inla.nonconvex.hull(coordinates(boundary.loc), 30),  # inner boundary
  inla.nonconvex.hull(coordinates(boundary.loc), 70)) # outer boundary

# Build the SPDE mesh
mesh <- inla.mesh.2d(boundary = bound, 
                     cutoff = 5, 
                     max.edge = c(30, 200),
                     min.angle=c(30, 21))

# Create SPDE object
spde <- inla.spde2.pcmatern(mesh, 
                            prior.range = c(50, 0.5),   # p(range < 50 = 0.5)
                            prior.sigma = c(0.1, 0.1))  # p(sigma > 0.1 = 0.1)


## Spatial field and data bookkeeping ----------

# Define data "groups" (in this case, years)
groups <- dat_stand$damage_year - (min(dat_stand$damage_year - 1))
ngroups <- length(unique(groups))

# Create projector matrix
# Use 'group' argument to tell INLA which observation
# belongs to which year
A <- inla.spde.make.A(mesh, 
                      group = groups,
                      loc = cbind(dat_stand$X_km, dat_stand$Y_km))

# Define the spatial fields
# This tells INLA which rows belong to the same year
field.z.idx <- inla.spde.make.index(name = 'z',
                                    n.spde = mesh$n,
                                    n.group = ngroups)
field.zc.idx <- inla.spde.make.index(name = 'zc',
                                     n.spde = mesh$n,
                                     n.group = ngroups)
field.u.idx <- inla.spde.make.index(name = 'u',
                                    n.spde = mesh$n,
                                    n.group = ngroups)


# Define data frames containing covariates 
N <- nrow(dat_stand)
X.z <- data.frame(z.Intercept = rep(1, N),
                  z.damage_year = dat_stand$damage_year,
                  z.Host_rich = dat_stand$Host_rich,
                  z.Host_BA = dat_stand$Host_BA,
                  z.Host_presence = dat_stand$Host_presence,
                  z.TWI = dat_stand$TWI,
                  z.HLI = dat_stand$HLI,
                  z.AET_norm = dat_stand$AET_norm,
                  z.tmin_norm = dat_stand$tmin_norm,
                  z.vpdmax_norm = dat_stand$vpdmax_norm,
                  z.tmin_anom = dat_stand$tmin_anom,
                  z.vpdmax_anom = dat_stand$vpdmax_anom)
X.y <- data.frame(y.Intercept = rep(1, N),
                  y.damage_year = dat_stand$damage_year,
                  y.Host_rich = dat_stand$Host_rich,
                  y.Host_BA = dat_stand$Host_BA,
                  y.Host_presence = dat_stand$Host_presence,
                  y.TWI = dat_stand$TWI,
                  y.HLI = dat_stand$HLI,
                  y.AET_norm = dat_stand$AET_norm,
                  y.tmin_norm = dat_stand$tmin_norm,
                  y.vpdmax_norm = dat_stand$vpdmax_norm,
                  y.tmin_anom = dat_stand$tmin_anom,
                  y.vpdmax_anom = dat_stand$vpdmax_anom)

# Make the data stacks for INLA model fitting
stack.z <- inla.stack(tag = "zobs",
                      data = list(y = cbind(as.vector(dat_stand$hotspot_bin), NA), link = 1),
                      A = list(A, 1),
                      effects = list(field.z.idx, X.z))
stack.y <- inla.stack(tag = "yobs",
                      data = list(y = cbind(NA, as.vector(dat_stand$hotspot_gt0)), link = 2),
                      A = list(A, 1),
                      effects = list(c(field.zc.idx, field.u.idx), X.y))
stack.all <- inla.stack(stack.z, stack.y)

# Save data stack as R object
saveRDS(stack.all, paste0("Objects/", prefix, "_stack.rds"))



## Model fitting ---------

# Slightly different model forms for each region

if(prefix == "SR"){
  
  form <- y ~ -1 + 
    f(z.damage_year, model = "rw1") +           # RW1 intercept
    f(y.damage_year, model = "rw1") +           # RW1 intercept
    z.Host_presence + z.Host_rich + z.Host_BA + # fixed effect of host cover
    y.Host_rich + y.Host_BA +                   # fixed effect of host cover
    z.TWI + z.HLI +                             # fixed effect of topography
    y.TWI + y.HLI +                             # fixed effect of topography
    z.AET_norm +                                # fixed effect of AET
    y.AET_norm +                                # fixed effect of AET
    z.tmin_norm * z.tmin_anom +                 # climate normals * weather anomalies
    z.vpdmax_norm * z.vpdmax_anom +             # climate normals * weather anomalies
    y.tmin_norm * y.tmin_anom +                 # climate normals * weather anomalies
    y.vpdmax_norm * y.vpdmax_anom +             # climate normals * weather anomalies
    f(z, model = spde, group = z.group,         # AR1 spatio-temporal random field for presence/absence
      control.group = list(model = "ar1")) +
    f(u, model = spde, group = u.group, 
      control.group = list(model = "ar1"))      # AR1 spatio-temporal random field for prevalence
}

if(prefix == "MR"){
  
  form <- y ~ -1 + 
    f(z.damage_year, model = "rw1") +           # RW1 intercept
    f(y.damage_year, model = "rw1") +           # RW1 intercept
    z.Host_presence + z.Host_rich + z.Host_BA + # fixed effect of host cover
    y.Host_rich + y.Host_BA +                   # fixed effect of host cover
    z.TWI + z.HLI +                             # fixed effect of topography
    y.TWI + y.HLI +                             # fixed effect of topography
    z.tmin_norm + z.vpdmax_norm + z.AET_norm +  # fixed effect of climate normals
    y.tmin_norm + y.vpdmax_norm + y.AET_norm +  # fixed effect of climate normals
    z.tmin_anom + z.vpdmax_anom +               # fixed effect of weather anomalies
    y.tmin_anom + y.vpdmax_anom +               # fixed effect of weather anomalies
    f(z, model = spde, group = z.group,         # AR1 spatio-temporal random field for presence/absence
      control.group = list(model = "ar1")) +
    f(u, model = spde, group = u.group, 
      control.group = list(model = "ar1"))      # AR1 spatio-temporal random field for prevalence
}

if(prefix == "C"){
  
  form <- y ~ -1 + 
    f(z.damage_year, model = "rw1") +           # RW1 intercept
    f(y.damage_year, model = "rw1") +           # RW1 intercept
    z.Host_presence + z.Host_rich + z.Host_BA + # fixed effect of host cover
    y.Host_rich + y.Host_BA +                   # fixed effect of host cover
    z.TWI + z.HLI +                             # fixed effect of topography
    y.TWI + y.HLI +                             # fixed effect of topography
    z.AET_norm + z.vpdmax_norm +                # fixed effect of climate normals
    y.AET_norm + y.vpdmax_norm +                # fixed effect of climate normals
    z.tmin_anom + z.vpdmax_anom +               # fixed effect of weather anomalies
    y.tmin_anom + y.vpdmax_anom +               # fixed effect of weather anomalies
    f(z, model = spde, group = z.group,         # AR1 spatio-temporal random field for presence/absence
      control.group = list(model = "ar1")) +
    f(u, model = spde, group = u.group, 
      control.group = list(model = "ar1"))      # AR1 spatio-temporal random field for prevalence
}

# Fit the model
# Zero-inflated binomial (type 0; hurdle model)
fit <- inla(form,
            family = c("binomial", "binomial"),
            Ntrials = c( rep(1, N), dat_stand$Ntrials_gt0 ),
            data = inla.stack.data(stack.all),
            control.compute = list(cpo = TRUE, dic  = TRUE, config = TRUE),
            control.predictor = list(A = inla.stack.A(stack.all), link = link),
            control.inla = list(strategy = "adaptive", int.strategy = "eb")) # for fast approximation

# Save model fit as R object
saveRDS(fit, paste0("Objects/", prefix, "_fit.rds"))

