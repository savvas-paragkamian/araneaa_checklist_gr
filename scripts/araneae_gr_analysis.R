#!/usr/bin/Rscript

## Script name: araneae_gr_analysis.R
##
## Purpose of script: 
## analysis of the occurrences data of the spiders of Greece
##
## Author: Savvas Paragkamian
##
## Date Created: 2024-11-05
##
library(tidyverse)
library(readxl)
library(sf)
library(terra)
library(units)
#library(vegan)

# load data 
hellenic_borders_shp <- sf::st_read("../data/hellenic_borders/hellenic_borders.shp")

## load data from excel file

araneae_gr_tax <- readxl::read_excel("../data/0_SpidoNETGR_02082024_Final.xlsx", sheet="Taxon Diversity")
araneae_gr_occ <- readxl::read_excel("../data/0_SpidoNETGR_02082024_Final.xlsx", sheet="Species Occurences")

names(araneae_gr_occ)[1] <- "Species"
names(araneae_gr_occ)[4] <- "decimalLongitude"
names(araneae_gr_occ)[3] <- "decimalLatitude"

### missing
which(is.na(araneae_gr_occ$decimalLongitude)) 

no_coordinates <- araneae_gr_occ[which(is.na(araneae_gr_occ$decimalLongitude)),]

write_delim(no_coordinates, "../results/no_coordinates.tsv",delim="\t")

### clean
### Agyneta pseudorurestris is twice.
#araneae_gr_tax |> group_by(Species) |> count(Species) |> arrange(desc(n)) 

### which species

#anti_species_occ <- araneae_gr_occ |> 
#    group_by(Species) |> 
#    summarise(occurrences=n()) |>
#    anti_join(araneae_gr_tax)

#write_delim(anti_species_occ, "../results/anti_species_occ.tsv",delim="\t")

#anti_species_tax <- araneae_gr_tax |>
#    anti_join(araneae_gr_occ)
    
#write_delim(anti_species_tax, "../results/anti_species_tax.tsv",delim="\t")

## there are 317 occurrences without coordinates

araneae_gr_occ_sf <- araneae_gr_occ |> 
    filter(!(is.na(decimalLongitude)), !(decimalLongitude>50)) |>
    st_as_sf(coords=c("decimalLatitude","decimalLongitude"),
             remove=F,
             crs="WGS84")


araneae_gr_inside_gr <- st_intersection(araneae_gr_occ_sf, hellenic_borders_shp)

## map

araneae_gr_base <- ggplot() +
    geom_sf(hellenic_borders_shp, mapping=aes()) +
    geom_point(araneae_gr_occ_sf,
            mapping=aes(x=decimalLongitude, y=decimalLatitude),
            color="orange",
            size=1.8,
            alpha=0.8,
            show.legend=T) +
    coord_sf(crs="WGS84") +
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "right",
          legend.box.background = element_blank())

ggsave("../figures/araneae_gr_base.png", 
       plot=araneae_gr_base, 
       height = 10, 
       width = 20,
       dpi = 600, 
       units="cm",
       device="png")

## taxonomic summary

araneae_gr_occ_tax <- araneae_gr_occ_sf |>
    left_join(araneae_gr_tax)

araneae_gr_occ_tax_f <- araneae_gr_occ_tax |> 
    st_drop_geometry() |>
    group_by(Species,Family) |>
    summarise(occurrences=n(), .groups="keep") |> 
    group_by(Family) |>
    summarise(n_species=n(),occurrences=sum(occurrences))

araneae_gr_occ_tax_f_l <- araneae_gr_occ_tax_f |>
    pivot_longer(-Family, names_to="Variable", values_to="number")

species_family_plot <- ggplot()+
  geom_col(data = araneae_gr_occ_tax_f_l,
           aes(x=Family, y= number, fill=Variable),
           width=0.82,
           position = position_dodge(width = 0.82),
           show.legend = T)+
  geom_text(data = araneae_gr_occ_tax_f_l,
            aes(x = Family ,y= number, label=number,group=Variable),
            position=position_dodge(width = 0.87),
            vjust=-0.25,
            size=5)+
#  scale_y_continuous(breaks = seq(0,280,20),
#                     limits = c(0,285),
#                     expand = c(0.01,0.4))+
#  scale_x_discrete(expand = c(0.01,0.4))+
#  scale_fill_manual(label=c("Caves","All species","Species endemic to Greece","Troglobiont species"),
#                    values = c("coral1","lightgoldenrod2","lightpink1","lightblue1"),
#                    name="")+
  labs(x="Family", y= "Count")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 18),
        axis.text.y=element_text(margin = margin(t = 0, r = 0, b = 0, l = 10,unit = "pt"),size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 18),
        axis.title = element_text(size=22),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 0.3),
        axis.line.y = element_line(colour = 'black', linewidth = 0.3),
        legend.position = c(0.13,0.87),
        legend.key.size = unit(1, "cm"))

ggsave("species_family_plot.png",
       plot = species_family_plot,
       device = "png",
       width = 60,
       height = 30,
       units = "cm",
       dpi = 100,
       path = "../figures/")


