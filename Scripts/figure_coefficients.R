# Figure 4: coefficients from INLA models

library(INLA)
library(tidyverse)

# Load model fit objects
fit_MR <- readRDS("Objects/MR_fit.rds")
fit_SR <- readRDS("Objects/SR_fit.rds")
fit_C <- readRDS("Objects/C_fit.rds")


## Summarize and combine fixed effects ------

# SR
fixed_SR <- summary(fit_SR)$fixed %>% 
  as.data.frame() %>% 
  mutate(Region = "Southern Rockies") %>%
  mutate(Important = ifelse(`0.025quant` * `0.975quant` < 0, "No", "Yes")) %>%
  mutate(parameter = rownames(.)) %>%
  mutate(response = str_sub(parameter, start = 1, end = 1)) %>%
  mutate(covariate = str_sub(parameter, start = 3, end = -1)) %>%
  mutate(covariate = str_remove(covariate, "z.")) %>%
  mutate(covariate = str_remove(covariate, "y.")) %>%
  mutate(response = ifelse(response == "z", "Occurrence", "Prevalence")) %>%
  mutate(response = factor(response, levels = c("Prevalence", "Occurrence"))) %>%
  mutate(covariate = factor(covariate, levels = c("Host_presence", "Host_BA", "Host_rich", 
                                                  "HLI", "TWI", "AET_norm",
                                                  "tmin_norm", "vpdmax_norm", 
                                                  "tmin_anom", "tmin_norm:tmin_anom",
                                                  "vpdmax_anom", "vpdmax_norm:vpdmax_anom")))

# MR
fixed_MR <- summary(fit_MR)$fixed %>% 
  as.data.frame() %>% 
  mutate(Region = "Middle Rockies") %>%
  mutate(Important = ifelse(`0.025quant` * `0.975quant` < 0, "No", "Yes")) %>%
  mutate(parameter = rownames(.)) %>%
  mutate(response = str_sub(parameter, start = 1, end = 1)) %>%
  mutate(covariate = str_sub(parameter, start = 3, end = -1)) %>%
  mutate(covariate = str_remove(covariate, "z.")) %>%
  mutate(covariate = str_remove(covariate, "y.")) %>%
  mutate(response = ifelse(response == "z", "Occurrence", "Prevalence")) %>%
  mutate(response = factor(response, levels = c("Prevalence", "Occurrence"))) %>%
  mutate(covariate = factor(covariate, levels = c("Host_presence", "Host_BA", "Host_rich", 
                                                  "HLI", "TWI", "AET_norm",
                                                  "tmin_norm", "vpdmax_norm", 
                                                  "tmin_anom", "tmin_norm:tmin_anom",
                                                  "vpdmax_anom", "vpdmax_norm:vpdmax_anom")))

# C
fixed_C <- summary(fit_C)$fixed %>% 
  as.data.frame() %>% 
  mutate(Region = "Cascades") %>%
  mutate(Important = ifelse(`0.025quant` * `0.975quant` < 0, "No", "Yes")) %>%
  mutate(parameter = rownames(.)) %>%
  mutate(response = str_sub(parameter, start = 1, end = 1)) %>%
  mutate(covariate = str_sub(parameter, start = 3, end = -1)) %>%
  mutate(covariate = str_remove(covariate, "z.")) %>%
  mutate(covariate = str_remove(covariate, "y.")) %>%
  mutate(response = ifelse(response == "z", "Occurrence", "Prevalence")) %>%
  mutate(response = factor(response, levels = c("Prevalence", "Occurrence"))) %>%
  mutate(covariate = factor(covariate, levels = c("Host_presence", "Host_BA", "Host_rich", 
                                                  "HLI", "TWI", "AET_norm",
                                                  "tmin_norm", "vpdmax_norm", 
                                                  "tmin_anom", "tmin_norm:tmin_anom",
                                                  "vpdmax_anom", "vpdmax_norm:vpdmax_anom")))


# All combined
fixed_all <- bind_rows(fixed_SR, fixed_MR, fixed_C) %>%
  mutate(Region = factor(Region, levels = c("Southern Rockies", "Middle Rockies", "Cascades")))

# Covariate labels
clabs <- tibble(covariate = c("Host_presence", "Host_BA", "Host_rich", 
                              "HLI", "TWI", "AET_norm",
                              "tmin_norm", "vpdmax_norm", 
                              "tmin_anom", "tmin_norm:tmin_anom",
                              "vpdmax_anom", "vpdmax_norm:vpdmax_anom"),
                covar_label = c("Host co-occurrence", "Host basal area", "Host richness",
                                "Topographic heat load index", "Topographic wetness index", "Annual AET (normal)",
                                "Winter min temp (normal)", "Summer max VPD (normal)",
                                "Winter min temp (anomaly)", "Winter min temp (norm:anom)",
                                "Summer max VPD (anomaly)", "Summer max VPD (norm:anom)"))
fixed_all <- left_join(fixed_all, clabs) %>%
  mutate(covar_label = factor(covar_label, levels = c("Host co-occurrence", "Host basal area", "Host richness",
                                                      "Topographic heat load index", "Topographic wetness index", "Annual AET (normal)",
                                                      "Winter min temp (normal)", "Summer max VPD (normal)",
                                                      "Winter min temp (anomaly)", "Winter min temp (norm:anom)",
                                                      "Summer max VPD (anomaly)", "Summer max VPD (norm:anom)")))

# Covariates that are 'inciting' factors
incite <- c('tmin_anom', 'vpdmax_anom', 
             'tmin_norm:tmin_anom', 'vpdmax_norm:vpdmax_anom')
fixed_all <- fixed_all %>%
  mutate(Predis_Incit = ifelse(covariate %in% incite, "Inciting factors", "Predisposing factors")) %>%
  mutate(Predis_Incit = factor(Predis_Incit, levels = c("Predisposing factors", "Inciting factors")))



# Plot ------

ggplot(fixed_all, aes(mean, covar_label, color = response, shape = Important)) +
  facet_grid(Predis_Incit~Region, scales = "free", space = "free_y") +
  geom_pointrange(aes(xmin = `0.025quant`, xmax = `0.975quant`),
                  position = position_dodge2(width = 0.4), size = 0.3, stroke = 0.5) +
  labs(x = "Standardized coefficient\n \n ", y = NULL) +
  geom_vline(xintercept = 0, linetype = "dotted", size = 0.5) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(name = NULL,
                     values = c(viridis::magma(1, begin = 0.8), "gray40")) +
  scale_shape_manual(values = c(1, 19)) +
  guides(color = guide_legend(reverse = TRUE,
                              override.aes = list(linetype = 0))) +
  guides(shape = "none") +
  theme_bw() + 
  theme(legend.position = c(0.4,-0.21),
        legend.direction = "horizontal",
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, vjust = 0),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 8),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray97"),
        plot.margin = margin(t = 5.5, l = 5.5, r = 5.5,
                             b = 25))

ggsave(paste0("Figures/coefficients.png"),
       width = 6.5, height = 4, units = "in", dpi = 350)

