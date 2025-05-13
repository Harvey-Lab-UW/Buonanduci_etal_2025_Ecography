# Figure 5: contribution of individual agents to hotspots

library(tidyverse)

# Read in data
agents <- read_csv("Data/agents.csv")
agents_contrib <- read_csv("Data/agents_hotspots_contrib.csv") %>%
  mutate(region = factor(region, levels = c("SR", "MR", "C")))

region_names = c(`SR` = "Southern\nRockies",
                 `MR` = "Middle\nRockies",
                 `C` = "Cascades")

# Plot
ggplot(filter(agents_contrib, code != 'hotspot'), 
       aes(y = code, x = contrib_prop)) +
  facet_grid(~region, labeller = as_labeller(region_names)) +
  geom_bar(stat = "identity") +
  scale_y_discrete(limits = rev, name = NULL,
                   labels = rev(str_to_sentence(agents$common_name))) +
  scale_x_continuous(name = "Proportion of hotspots") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 8, vjust = 0),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 8),
        panel.grid = element_blank())

ggsave("Figures/agents_hotspots_contrib.png",
       width = 4, height = 3, units = "in", dpi = 350)

