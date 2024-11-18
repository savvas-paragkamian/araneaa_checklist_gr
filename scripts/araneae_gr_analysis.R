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
municipalities_names <- readxl::read_excel("../data/municipalities_shape_file/names_municipalities_gr_eng.xlsx") |>
    mutate(KWD_YPES=as.character(KWD_YPES))

hellenic_municipalities_shp <- sf::st_read("../data/municipalities_shape_file/municipalities_Kallikratis_plan_Greece.shp") |>
    left_join(municipalities_names) |>
    dplyr::select(-NAME)

# https://gadm.org
greece_regions <- sf::st_read("../data/gadm41_GRC_shp/gadm41_GRC_2.shp")

#greece_regions <- c("Athos","East Macedonia and Thrace","Attica ","West Greece","West Macedonia","Ionian Islands ","Epirus ","Central Macedonia","Crete","South Aegean","Peloponnese ","Central Greece ","Thessaly","North Aegean")


## load data from excel file

araneae_gr_tax <- readxl::read_excel("../data/spiders_gr_dataset.xlsx", sheet="Taxon")
araneae_gr_occ <- readxl::read_excel("../data/spiders_gr_dataset.xlsx", sheet="Occurrence")
araneae_gr_ref <- readxl::read_excel("../data/spiders_gr_dataset.xlsx", sheet="References")


### missing
which(is.na(araneae_gr_occ$decimalLongitude)) 

no_coordinates <- araneae_gr_occ[which(is.na(araneae_gr_occ$decimalLongitude)),]

write_delim(no_coordinates, "../results/no_coordinates.tsv",delim="\t")

### clean
### Agyneta pseudorurestris is twice.
#araneae_gr_tax |> group_by(scientificName) |> count(scientificName) |> arrange(desc(n)) 

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

# keep only the occurrences inside Greece. 85 occurrences are removed
# this takes 3-5 minutes to run
araneae_gr <- st_intersection(araneae_gr_occ_sf,greece_regions)

## map

araneae_gr_base <- ggplot() +
    geom_sf(greece_regions, mapping=aes(),color="gray80") +
    geom_point(araneae_gr,
            mapping=aes(x=decimalLatitude, y=decimalLongitude),
            color="orange",
            size=1,
            alpha=0.5,
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
       height = 20, 
       width = 20,
       dpi = 300, 
       units="cm",
       device="png")


## taxonomic summary

araneae_gr_occ_tax <- araneae_gr |>
    left_join(araneae_gr_tax,
              by=c("scientificName"="scientificName"))


## IUCN categories per family

araneae_gr_occ_tax_iucn <- araneae_gr_occ_tax |>
    distinct(family,scientificName,iucnStatus) |>
    group_by(family,iucnStatus) |> 
    summarise(species=n(), .groups="keep") |> 
    ungroup() |>
    pivot_wider(names_from=iucnStatus,
                values_from=species,
                values_fill=0)


## endemics per family
araneae_endemics_family <- araneae_gr_occ_tax |>
    distinct(family,scientificName,endemic) |>
    na.omit(endemic) |>
    group_by(family) |>
    summarise(endemics=n())

## summary per family
###     distinct(family, scientificName,decimalLongitude,decimalLatitude) 
###     because some references have the same occurrence

araneae_gr_occ_tax_f <- araneae_gr_occ_tax |> 
    st_drop_geometry() |>
    distinct(family, scientificName,decimalLongitude,decimalLatitude) |>
    group_by(scientificName,family) |>
    summarise(occurrences=n(), .groups="keep") |> 
    group_by(family) |>
    summarise(species=n(),occurrences=sum(occurrences)) |>
    left_join(araneae_endemics_family) |>
    mutate(across(everything(), ~ replace_na(.x, 0)))

araneae_f <- araneae_gr_occ_tax_f |>
    left_join(araneae_gr_occ_tax_iucn)

araneae_family_summary <- araneae_f |>
    mutate(family="total") |>
    group_by(family) |>
    summarise(across(everything(), ~ sum(.x,na.rm = TRUE))) |>
    bind_rows(araneae_f)

write_delim(araneae_family_summary,
            "../results/araneae_family_summary.tsv",
            delim="\t")

araneae_gr_occ_tax_f_l <- araneae_gr_occ_tax_f |>
    pivot_longer(-family, names_to="Variable", values_to="number")

## bar plot of families

species_family_plot <- ggplot()+
  geom_col(data = araneae_gr_occ_tax_f_l,
           aes(x=family, y= number, fill=Variable),
           width=0.82,
           position = position_dodge(width = 0.82),
           show.legend = T)+
  geom_text(data = araneae_gr_occ_tax_f_l,
            aes(x = family ,y= number, label=number,group=Variable),
            position=position_dodge(width = 0.87),
            vjust=-0.25,
            size=3)+
  scale_y_continuous(breaks = seq(0,4000,500),
                     limits = c(0,4000))+
  scale_x_discrete(expand = c(0.01,0.4))+
  scale_fill_manual(values = c("species"="darkseagreen4",
                               "occurrences"="darkorange",
                               "endemics"="firebrick1"),
                    name="")+
  labs(x="Family", y= "Count")+
  guides(fill= guide_legend(position = "inside"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        axis.text.y=element_text(margin = margin(t = 0, r = 0, b = 0, l = 10,unit = "pt"),size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.title = element_text(size=16),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = 'black', linewidth = 0.3),
        axis.line.y = element_line(colour = 'black', linewidth = 0.3),
        legend.position.inside = c(0.08,0.82),
        legend.key.size = unit(1, "cm"))

ggsave("species_family_plot.png",
       plot = species_family_plot,
       device = "png",
       width = 60,
       height = 20,
       units = "cm",
       dpi = 300,
       path = "../figures/")

## summary per region
# occurrences
region_summary_occ <- araneae_gr_occ_tax |> 
    distinct(NAME_2,family, scientificName,endemic,decimalLongitude,decimalLatitude) |>
    group_by(NAME_2) |>
    summarise(occurrences=n())

# species
region_summary_sp <- araneae_gr_occ_tax |> 
    distinct(NAME_2,scientificName) |>
    group_by(NAME_2) |>
    summarise(species=n())

# endemics
region_summary_end <- araneae_gr_occ_tax |> 
    distinct(NAME_2,scientificName, endemic) |>
    na.omit(endemic) |> 
    group_by(NAME_2) |>
    summarise(endemics=n())

region_all <- region_summary_occ |>
    left_join(region_summary_sp) |>
    left_join(region_summary_end)

write_delim(region_all, "../results/regions_all_stats.tsv", delim="\t")

# References summary

araneae_gr_occ_ref <- araneae_gr_occ_tax |>
    left_join(araneae_gr_ref,
              by=c("associatedReferences"="associatedReferences"))

araneae_ref_year <- araneae_gr_occ_ref |> 
    distinct(associatedReferences,year) |>
    group_by(year) |>
    summarise(occurrance=n()) |>
    arrange(year) |>
    mutate(Cumulative_occurrance=cumsum(occurrance),
           Classification="Publications") 


## Species knowledge accumulation

### all species per year per reference
species_cumulative <- araneae_gr_occ_ref |> 
    distinct(scientificName, associatedReferences,year) |>
    arrange(year) |>
    mutate(Duplicates=duplicated(scientificName)) |> 
    mutate(First_occurrance=if_else(Duplicates=="FALSE",1,0)) |> 
    na.omit() |> 
    filter(First_occurrance==1) |>
    group_by(year) |> 
    summarise(species_year= n()) |>
    arrange(year) |> 
    mutate(Cumulative_occurrance= cumsum(species_year)) |> 
    mutate(Classification="All species") |> 
    dplyr::select(-c(species_year)) |>
    distinct()

### endemic species accumulation across the years
endemic_cumulative_species <- araneae_gr_occ_ref |> 
    distinct(scientificName,endemic, associatedReferences,year) |>
    arrange(year) |>
    mutate(Duplicates=duplicated(scientificName)) |> 
    mutate(First_occurrance=if_else(Duplicates=="FALSE",1,0)) |> 
    na.omit() |> 
    filter(First_occurrance==1) |>
    group_by(year,endemic) |> 
    summarise(endemic_species_year= n(), .groups="keep") |>
    ungroup() |>
    arrange(year) |> 
    mutate(Cumulative_occurrance= cumsum(endemic_species_year)) |>
    mutate(Classification="Endemic species to Greece") |> 
    dplyr::select(-c(endemic,endemic_species_year)) |> 
    ungroup() 

### combine data
araneae_accumulation <- bind_rows(araneae_ref_year, species_cumulative,endemic_cumulative_species)

### timeline figure
araneae_accumulation_plot <- ggplot()+
    geom_line(data=araneae_accumulation,
              aes(x=year,
                  y= Cumulative_occurrance,color=Classification),
              linewidth=1,
              show.legend = T)+
    #ggtitle("Class")+
    scale_x_continuous(breaks = seq(1830,2030,10),
                       limits = c(1830,2030),
                       expand=c(0.015,0))+
    scale_y_continuous(breaks = seq(0,1400,100),
                       limits = c(0,1400),
                       expand = c(0.01,0))+
    scale_color_manual(values =c("Endemic species to Greece"="firebrick1",
                                 "All species"="darkseagreen4"))+
    labs(x="Years",
         y="Cumulative number of species")+
    theme_bw()+
    theme(panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          legend.text = element_text(size = 18),
          axis.text.y=element_text(margin = margin(t = 0, r = 5, b = 0, l = 15,unit = "pt"),size = 18),
          axis.text.x = element_text(margin = margin(t = 5, r = 0, b = 15, l = 0,unit = "pt"),size = 18),
          axis.title = element_text(size=22),
          panel.border = element_blank(),
          axis.line.x = element_line(colour = 'black', linewidth = 0.3), 
          axis.line.y = element_line(colour = 'black', linewidth = 0.3),
          legend.position = c(0.13,0.87), 
          legend.key.size = unit(1.5, "cm"), 
          legend.title = element_blank())
  
ggsave("araneae_accumulation_plot.png",
       plot = araneae_accumulation_plot,
       device = "png",
       width = 50,
       height = 36,
       units = "cm",
       dpi = 300,
       path = "../figures/")

