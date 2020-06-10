
library(ggplot2)
library(bcmaps)
library(bcmapsdata)
library(tidyverse)

#Get Eco_Prov as geometry, exclude 'NEP'
ecoprov <- ecoprovinces() %>% 
  dplyr::filter(ECOPROVINCE_CODE!='NEP')


#Get boundry of processing extent
pro_extnt<- read_sf("Manuscript/tatolatex/Figures/Bev_SD_Extnt.gpkg")


#A rabbit hole 
# bcn <- bc_neighbours %>% mutate(name = name, type = "neighbours") %>% select(name, type)
# pro <- pro_extnt %>% mutate(name = "extent", type = "extent") %>% select(name, type)
# eco <- ecoprovinces %>% mutate(name = ECOPROVINCE_NAME, type = "ecoprov") %>% select(name, type)
# 
# rbind(bcn, pro)


ggplot() +
    geom_sf(data = bc_neighbours, aes(fill = name)) +
    geom_sf(data = pro_extnt, fill = NA, color = "black", size = 3, linetype = 2) +
    geom_sf(data = ecoprov, fill = NA) +
    geom_sf_label(data = ecoprovinces, aes(label = ECOPROVINCE_NAME)) +
    labs(x = "", y = "", fill = "Regions")


ggsave(filename = "Manuscript/tatolatex/Figures/study_area.png", device = 'png')



library(RColorBrewer)
library(bcmaps)
library(sf)
library(boot)
library(gridExtra)
library(raster)
library(RStoolbox)
library(ggrepel)
library(ggridges)
library(scales)
library(egg)

dem_msk <- raster("../snow_summaries_note/Data/hillshade_clip_250m.tif")

ggR(dem_msk)


#Make plot of NRDs, save to working dir

  # TRIM Hillshade (optional)
ggplot() + 
  geom_sf(data = bc_neighbours, aes(fill = iso_a2), show.legend = F) + 
  scale_fill_manual(values = c("grey","#AADAFF","dark grey")) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  ggR(dem_msk, ggLayer = T, maxpixels = ncell(dem_msk)/4) +
  geom_sf(data = ecoprov, color = "grey30", size = 0.5, fill = NA) +
  # NR District names to legend (without outline)
  geom_sf(data = pro_extnt, fill = NA, size = 2 ,linetype = 2, color = "white") +
  geom_label_repel(data = ecoprov,
                   aes(label = sub(pattern = " ", replacement = "\n", str_to_title(ECOPROVINCE_NAME)), geometry = geometry),
                   fontface = "bold",
                   stat = "sf_coordinates",
                   point.padding = NA,
                   arrow = NULL, #arrow(length = unit(0.02, "npc")),
                   segment.size = 0.5, segment.colour = NA,
                   nudge_x = 0,
                   nudge_y = 0,
                   hjust = 0.5,
                   vjust = 0.5,
                   show.legend = F) +
  geom_sf_text(data = filter(bc_neighbours, postal %in% c("BC","AB","YT","AK", "WA")), aes(label = str_to_upper(name)), show.legend = F, fontface = "bold", size = 6) + 
  geom_sf_text(data = filter(bc_neighbours, postal %in% c("BC","AB","YT","AK", "WA")), aes(label = str_to_upper(name)), show.legend = F, fontface = "bold", 
               size = 6, 
               color = "white",
               nudge_x = 3000,
               nudge_y = 3000) + 
  theme_minimal() +
  labs(x = "", y = "")

ggsave(filename = "Manuscript/tatolatex/Figures/study_area.png", device = 'png', width = 9, height = 7.5)


