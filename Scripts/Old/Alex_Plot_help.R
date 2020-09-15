
library(tidyverse)

df <- lda_can$coeffs.raw

df %>% 
  as.data.frame() %>% 
  mutate(names = row.names(lda_can$coeffs.raw)) %>% 
  ggplot() + 
    geom_point(aes(Can1, Can2), size = 3) +
    geom_text(aes(Can1, Can2, label = names), nudge_x = 0.005, nudge_y = 0.005) + 
    theme_bw() + 
    theme(aspect.ratio = 1) +
    geom_abline(slope = 1, intercept = 0) +
    geom_vline(xintercept = 0, linetype = 2) +
    labs(x = "This is x", y = " this is Y", title = "Title")

ggsave(filename = "test.pdf", width = 8, height = 8)
getwd()


library(bcmaps)
library(ggspatial)
ggplot() +
  # layer_spatial(data = pc1_ld) +
  geom_sf(data = bc_bound(), fill = NA, aes(color = "British Columbia")) +
  geom_sf(data = ecoprovinces(), fill = NA, aes(color = "Eco-Provinces")) +
  geom_sf_label(data = ecoprovinces(), aes(label = ECOPROVINCE_CODE, fill = ECOPROVINCE_CODE)) +
  scale_color_manual(values = c("black","dark grey")) + 
  scale_fill_discrete()

#Git test
#Another Change
