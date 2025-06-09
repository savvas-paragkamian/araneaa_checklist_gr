#!/usr/bin/env Rscript

## Script name: araneae_gr_analysis.R
##
## Purpose of script: 
## analysis of the occurrences data of the spiders of Greece
##
## Author: Savvas Paragkamian
##
## Date Created: 2024-11-05
##
library(dplyr)
library(readr)
library(ggplot2)
library(sf)

greece_regions <- sf::st_read("../data/gadm41_GRC_shp/gadm41_GRC_2.shp")

## opiliones

opi_data <- read_delim("../data/Iczegar_simplemaps_opi.csv",delim=",")

################################## quality check ###################################
sort(unique(opi_data$sci_name))

# load to sf object to check coordinates
points <- opi_data |> 
    st_as_sf(coords=c("dec_long","dec_lat"),
             remove=F,
             crs="WGS84")

# which are inside land and hellenic borders

dist_matrix <- st_distance(points, greece_regions)
points$min_dist_m <- apply(dist_matrix, 1, min)

# points inside
points_inside_or_touching <- points[points$min_dist_m ==0, ]

# Keep all points away from shore m of any polygon
points_away <- points[points$min_dist_m > 0 ,]

## save points outside
points_outside_cols <- points |>
    filter(min_dist_m >0) 

write_delim(points_outside_cols,
            "../results/opiliones__outside_gr_land.tsv",
            delim="\t")



################################## maps ###################################
## map all

opi_gr_base <- ggplot() +
    geom_sf(greece_regions, mapping=aes(),color="gray70") +
    geom_point(opi_data,
            mapping=aes(y=dec_lat, x=dec_long,color=Endemic),
            size=1,
            alpha=0.7,
            show.legend=T) +
    coord_sf(crs="WGS84") +
    ggtitle("Opiliones of Greece")+
    scale_color_manual(values = c("yes"="firebrick1",
                               "no"="cadetblue"))+
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "inside",
          legend.position.inside = c(0.9,0.9),
          legend.box.background = element_blank())

ggsave("../figures/opiliones_gr_base.png", 
       plot=opi_gr_base, 
       height = 20, 
       width = 20,
       dpi = 300, 
       units="cm",
       device="png")

## maps with points away

opi_gr_base_away <- ggplot() +
    geom_sf(greece_regions, mapping=aes(),color="gray70") +
    geom_point(points_outside_cols,
            mapping=aes(y=dec_lat, x=dec_long,color="Away"),
            size=1,
            alpha=0.7,
            show.legend=T) +
    geom_point(points_inside_or_touching,
            mapping=aes(y=dec_lat, x=dec_long,color="Inside"),
            size=1,
            alpha=0.7,
            show.legend=T) +
    coord_sf(crs="WGS84") +
    ggtitle("Opiliones of points away and inside Greece")+
    scale_color_manual(values = c("Away"="firebrick1",
                               "Inside"="cadetblue"))+
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "inside",
          legend.position.inside = c(0.9,0.9),
          legend.box.background = element_blank())

ggsave("../figures/opiliones_gr_aways.png", 
       plot=opi_gr_base_away, 
       height = 20, 
       width = 20,
       dpi = 300, 
       units="cm",
       device="png")
## opiliones redlist
## ICZEGAR conference in Athens
iucn_colors <- c(
  "DD" = "gray42",    # Data Deficient
  "LC" = "#00FF00",    # Least Concern
  "NT" = "#9ACD32",    # Near Threatened
  "VU" = "#FFFF00",    # Vulnerable
  "EN" = "#FF9900",    # Endangered
  "CR" = "#FF0000"    # Critically Endangered
)

opi_data$IUCN_cat2 <- factor(opi_data$IUCN_cat, levels = c("DD", "LC", "NT", "VU", "EN", "CR"))

opi_data_tax <- opi_data[order(opi_data$IUCN_cat2), ]

opi_gr_redlist <- ggplot() +
    geom_sf(greece_regions, mapping=aes(),color="gray70") +
    geom_point(opi_data_tax,
            mapping=aes(y=dec_lat, x=dec_long,color=IUCN_cat2),
            size=1.2,
            alpha=1,
            show.legend=T) +
    coord_sf(crs="WGS84") +
    scale_color_manual(values =iucn_colors, name="IUCN category")+
    theme_bw()+
    theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.title = element_text(size = 8),
        legend.position = "inside",
        legend.position.inside = c(0.9, 0.8),
        legend.box.background = element_blank()
  )

ggsave("../figures/opiliones_gr_redlist.png", 
       plot=opi_gr_redlist, 
       height = 20, 
       width = 26,
       dpi = 300, 
       units="cm",
       device="png")

### facet

opi_gr_redlist_f <- opi_gr_redlist +
    facet_wrap(~IUCN_cat2,nrow=2)+
    theme(
        legend.position = "inside",
        legend.position.inside = c(0.95, 0.85)
    )


ggsave("../figures/opiliones_gr_redlist_f.png", 
       plot=opi_gr_redlist_f, 
       height = 20, 
       width = 35,
       dpi = 300, 
       units="cm",
       device="png")


