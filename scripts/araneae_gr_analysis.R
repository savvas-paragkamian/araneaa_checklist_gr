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
library(ggpubr)
library(readxl)
library(sf)
library(terra)
library(units)
#library(vegan)

# load data 
gr_1km <- sf::st_read("../data/eea_1km/gr_1km.shp") |>
    st_transform(., crs="WGS84")

# https://gadm.org
greece_regions <- sf::st_read("../data/gadm41_GRC_shp/gadm41_GRC_2.shp")
greece_municipalities <- sf::st_read("../data/gadm41_GRC_shp/gadm41_GRC_3.shp")

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
# add the cellcode of the EEA reference grid
araneae_gr <- st_intersection(araneae_gr_occ_sf,greece_regions) |>
    st_join(gr_1km)

#araneae_gr_1km <- st_join(araneae_gr,gr_1km)

################################## islands #############################
###
# Explode multipolygons into individual polygons
gadm_single <- greece_municipalities |>
    st_make_valid() |>
    st_cast("POLYGON", do_split = TRUE) |>
    mutate(id = row_number())  # Add unique IDs to track original geometries
# Check connectivity of polygons

#connectivity <- st_relate(gadm_single, gadm_single, pattern = "****0****")
connectivity <- st_intersects(gadm_single, gadm_single)

# Determine islands (features with no neighbors are islands)
gadm_single$is_island <- sapply(connectivity, function(x) length(x) == 1)

greece_islands <- gadm_single |>
  mutate(is_island= ifelse(NAME_2 %in% c("North Aegean","Crete","South Aegean","Ionian Islands"),TRUE,is_island)) |>
  mutate(region_type = ifelse(is_island, "Island", "Mainland")) |>
  dplyr::select(is_island,region_type,NAME_2, NAME_3,id) |>
  mutate(area_island=round(set_units(st_area(geometry),km^2),4)) 

# Crete and Evia ara the only islands in greece that has multiple municipalities
# so the mainland of Crete is filtered to join all municipalities
# and keep satelite islands separate.
crete_only <- greece_islands |>
    filter(NAME_2=="Crete") |>
    arrange(desc(area_island)) |>
    filter(area_island>set_units(100,km^2))

crete_only_one <- st_union(crete_only) |>
    sf::st_as_sf() |>
    rename("geometry"="x") |>
    mutate(is_island=TRUE, region_type="Island", NAME_2="Crete", NAME_3="All Crete", id=0) |> 
    mutate(area_island=round(set_units(st_area(geometry),km^2),4)) 

# Evia
evia <- greece_regions |> filter(NAME_2=="Central Greece") |>
    st_make_valid() |>
    st_cast("POLYGON", do_split = TRUE) |>
    mutate(id = row_number()) |> 
    mutate(area_island=round(set_units(st_area(geometry),km^2),4)) |>
    filter(area_island>set_units(1000,km^2)) |>
    filter(area_island<set_units(4000,km^2)) |>
    dplyr::select(area_island) |>
    mutate(is_island=TRUE, region_type="Island", NAME_2="Central Greece", NAME_3="Evia", id=1000) 

# first the crete mainland polygons are excluded and then the whole crete is
# binded to the other
greece_islands_final <- greece_islands |>
    filter(!(id %in% crete_only$id)) |>
    bind_rows(crete_only_one) |>
    bind_rows(evia)

araneae_gr_all <- st_intersection(araneae_gr,greece_islands_final)

islands_gr <- ggplot() +
  geom_sf(greece_islands_final,mapping=aes(fill = region_type)) +
    #geom_point(araneae_gr_a,
    #        mapping=aes(x=decimalLatitude, y=decimalLongitude,color=is_island),
    #        size=1,
    #        alpha=0.7,
    #        show.legend=T) +
  theme_bw() +
  labs(title = "Islands vs Mainland in GADM Data")

ggsave("islands_gr.png",
       plot = islands_gr,
       device = "png",
       width = 20,
       height = 20,
       units = "cm",
       dpi = 300,
       path = "../figures/")

############ to do ##########
### which are the insular only species


############################ taxonomic summary ###########################
## master dataframe
araneae_gr_occ_tax <- araneae_gr_all |>
    left_join(araneae_gr_tax,
              by=c("scientificName"="scientificName")) |>
    mutate(taxonDistribution=if_else(is.na(endemic),"Non endemic","Endemic")) 

write_delim(araneae_gr_occ_tax,"../results/araneae_gr_occ_tax.tsv",delim="\t")

##  taxonomy on islands

island_occurrences <- araneae_gr_occ_tax |>
    distinct(scientificName,decimalLongitude,decimalLatitude,taxonDistribution,region_type) |>
    group_by(region_type,taxonDistribution) |>
    summarise(occurrences=n(),.groups="keep")

region_type_species <- araneae_gr_occ_tax |>
    distinct(scientificName,decimalLongitude,decimalLatitude,region_type,taxonDistribution) |>
    group_by(scientificName,region_type,taxonDistribution) |>
    summarise(occurrences=n(), .groups="keep")

singular_region_species <- region_type_species |>
    group_by(scientificName) |>
    filter(n() == 1) |>
    ungroup() 

duplicate_region_species <- region_type_species |>
    group_by(scientificName) |>
    filter(n() == 2) |>
    ungroup() 

singular_region_species_summary <- singular_region_species |>
    group_by(taxonDistribution,region_type) |>
    summarise(singular_species=n(), .groups="keep") |>
    left_join(island_occurrences, by=c("region_type","taxonDistribution"))

write_delim(singular_region_species_summary,"../results/singular_region_species_summary.tsv",delim="\t")

#### islands stats

island_stats <- araneae_gr_occ_tax |>
    distinct(scientificName,taxonDistribution,region_type, NAME_2, id,area_island) |>
    filter(region_type=="Island") |>
    group_by(id,area_island,NAME_2,region_type,taxonDistribution) |>
    summarise(species=n(), .groups="keep")

### island endemics plots
island_stats_plot <- ggplot()+
    geom_point(island_stats,
               mapping=aes(x=area_island,
                   y=species,
                   color=taxonDistribution))+

    theme_bw()

ggsave("island_stats_endemics.png",
       plot = island_stats_plot,
       device = "png",
       width = 30,
       height = 30,
       units = "cm",
       dpi = 300,
       path = "../figures/")

island_stats_log_plot <- ggplot()+
    geom_point(island_stats,
               mapping=aes(x=log10(area_island),
                   y=log10(species),
                   color=taxonDistribution))+

    theme_bw()

ggsave("island_stats_log_endemics.png",
       plot = island_stats_log_plot,
       device = "png",
       width = 30,
       height = 30,
       units = "cm",
       dpi = 300,
       path = "../figures/")

### island types plots
island_stats_pelagos_plot <- ggplot()+
    geom_point(island_stats,
               mapping=aes(x=area_island,
                   y=species,
                   color=NAME_2))+

    theme_bw()

ggsave("island_stats_pelagos.png",
       plot = island_stats_pelagos_plot,
       device = "png",
       width = 30,
       height = 30,
       units = "cm",
       dpi = 300,
       path = "../figures/")

island_stats_pelagos_log_plot <- ggplot()+
    geom_point(island_stats,
               mapping=aes(x=log10(area_island),
                   y=log10(species),
                   color=NAME_2))+
#    scale_color_brewer(palette = "Set2")+
    theme_bw()

ggsave("island_stats_pelagos_log.png",
       plot = island_stats_pelagos_log_plot,
       device = "png",
       width = 30,
       height = 30,
       units = "cm",
       dpi = 300,
       path = "../figures/")

island_stats_pelagos_log_plot_f <- island_stats_pelagos_log_plot + facet_wrap(~taxonDistribution)


ggsave("island_stats_pelagos_log_f.png",
       plot = island_stats_pelagos_log_plot_f,
       device = "png",
       width = 30,
       height = 15,
       units = "cm",
       dpi = 300,
       path = "../figures/")

# From the total 1299 species, the 735 appear only in mainland Greece or in greek islands.
#


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

araneae_endemics_all_all <- araneae_gr_occ_tax |> 
    distinct(scientificName,
             family,
             decimalLongitude,decimalLatitude,
             endemic,
             iucnStatus,CELLCODE) |>
    filter(!is.na(endemic)) 
## all endemics
araneae_endemics_all <- araneae_gr_occ_tax |> 
    distinct(scientificName,
             family,
             decimalLongitude,decimalLatitude,
             endemic,
             iucnStatus,CELLCODE) |>
    filter(!is.na(endemic)) |>
    group_by(scientificName)|>
    mutate(occurrences=n(),
           sq_km=length(unique(CELLCODE)),
           species=1) |>
    distinct(scientificName,species,occurrences,iucnStatus,sq_km) |>
    ungroup()

araneae_endemics_iucn <- araneae_endemics_all |>
    dplyr::select(iucnStatus,species) |>
    group_by(iucnStatus) |>
    summarise(species=sum(species)) |>
    pivot_wider(names_from=iucnStatus,
                values_from=species,
                values_fill=0)

araneae_endemics_sum <- araneae_endemics_all |>
    dplyr::select(species,occurrences,sq_km) |>
    mutate(family="endemics") |>
    group_by(family) |>
    summarise(across(everything(), ~ sum(.x,na.rm = TRUE))) |>
    bind_cols(araneae_endemics_iucn)

## eea sq km per family
##



## summary per family
###     distinct(family, scientificName,decimalLongitude,decimalLatitude) 
###     because some references have the same occurrence

araneae_gr_occ_tax_f <- araneae_gr_occ_tax |> 
    st_drop_geometry() |>
    distinct(family,
             scientificName,
             decimalLongitude,
             decimalLatitude, CELLCODE) |>
    group_by(scientificName,family) |>
    summarise(occurrences=n(),
              sq_km=length(unique(CELLCODE)),
              .groups="keep") |> 
    group_by(family) |>
    summarise(species=n(),
              occurrences=sum(occurrences),
              sq_km=sum(sq_km)) |>
    left_join(araneae_endemics_family) |>
    mutate(across(everything(), ~ replace_na(.x, 0)))

araneae_f <- araneae_gr_occ_tax_f |>
    left_join(araneae_gr_occ_tax_iucn)

## total

araneae_total_iucn <- araneae_f  |>
    dplyr::select(-c(occurrences,sq_km)) |>
    mutate(family="total") |>
    group_by(family) |>
    summarise(across(everything(), ~ sum(.x,na.rm = TRUE))) 

araneae_total_area <- araneae_gr_occ_tax |>
    mutate(family="total") |>
    distinct(decimalLongitude,decimalLatitude,CELLCODE, family) |>
    group_by(family) |>
    summarise(occurrences=n(), sq_km=length(unique(CELLCODE)))

araneae_family_summary <- araneae_total_area |>
    left_join(araneae_total_iucn) |>
    bind_rows(araneae_endemics_sum,araneae_f)

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
                               "sq_km"="cadetblue",
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

################################## gbif ###################################
gbif_arachnida <- read_delim("../data/0058835-241126133413365/occurrence.txt", delim="\t")

# filter
gbif_f <- gbif_arachnida |>
    dplyr::distinct(phylum,
                  class,
                  order,
                  family,
                  species,
                  acceptedScientificName,
                  collectionCode,
                  basisOfRecord,
                  institutionCode,
                  year,
                  decimalLatitude,
                  decimalLongitude)


gbif_occ_citizen <- gbif_f |>
    filter(!is.na(decimalLatitude)) |>
    filter(basisOfRecord=="HUMAN_OBSERVATION") |>
    st_as_sf(coords=c("decimalLongitude", "decimalLatitude"),
             remove=F,
             crs="WGS84") |>
    filter(order=="Araneae") |>
    st_intersection(greece_regions) |>
    st_join(gr_1km)


gbif_occ <- gbif_f |>
    filter(!is.na(decimalLatitude)) |>
    filter(basisOfRecord!="HUMAN_OBSERVATION") |>
    st_as_sf(coords=c("decimalLongitude", "decimalLatitude"),
             remove=F,
             crs="WGS84")

gbif_araneae_gr <- gbif_occ |>
    filter(order=="Araneae") |>
    st_intersection(greece_regions) |>
    st_join(gr_1km)


gbif_join <- bind_rows(gbif_occ_citizen,gbif_araneae_gr)
################################## maps ###################################
## map

araneae_gr_base <- ggplot() +
    geom_sf(greece_regions, mapping=aes(),color="gray70") +
    geom_point(araneae_gr_occ_tax,
            mapping=aes(x=decimalLatitude, y=decimalLongitude,color=taxonDistribution),
            size=1,
            alpha=0.7,
            show.legend=T) +
    coord_sf(crs="WGS84") +
    ggtitle("Araneae of Greece")+
    scale_color_manual(values = c("Endemic"="firebrick1",
                               "Non endemic"="cadetblue"))+
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "inside",
          legend.position.inside = c(0.9,0.9),
          legend.box.background = element_blank())

ggsave("../figures/araneae_gr_base.png", 
       plot=araneae_gr_base, 
       height = 20, 
       width = 20,
       dpi = 300, 
       units="cm",
       device="png")

## gbif
araneae_gr_gbif <- ggplot() +
    geom_sf(greece_regions, mapping=aes(),color="gray70") +
    geom_point(gbif_araneae_gr,
            mapping=aes(y=decimalLatitude, x=decimalLongitude,color=basisOfRecord),
            size=1,
            alpha=0.7,
            show.legend=T) +
    coord_sf(crs="WGS84") +
    ggtitle("GBIF Araneae of Greece")+
    scale_color_manual(values = c(
                                  #"HUMAN_OBSERVATION"="springgreen2",
                                  "PRESERVED_SPECIMEN"="lightsalmon2",
                                  "MATERIAL_CITATION"="darkorchid2",
                                  "OCCURRENCE"="springgreen4",
                                  "MATERIAL_SAMPLE"="khaki3"
                               ))+
    theme_bw()+
    theme(axis.title=element_blank(),
          axis.text=element_text(colour="black"),
          legend.title = element_text(size=8),
          legend.position = "inside",
          legend.position.inside = c(0.85,0.8),
          legend.box.background = element_blank())

ggsave("../figures/araneae_gr_gbif.png", 
       plot=araneae_gr_gbif, 
       height = 20, 
       width = 20,
       dpi = 300, 
       units="cm",
       device="png")



## families
family_names <- unique(araneae_gr_occ_tax$family)
family_plots <- list()

for (i in seq_along(family_names)) {

    print(family_names[i])

    family_data <- araneae_gr_occ_tax |>
        filter(family==family_names[i])

    araneae_family <- ggplot() +
        geom_sf(greece_regions, mapping=aes(),color="gray70") +
        geom_point(family_data,
                mapping=aes(x=decimalLatitude, y=decimalLongitude,color=taxonDistribution),
                size=1,
                alpha=0.7,
                show.legend=T) +
        coord_sf(crs="WGS84") +
        ggtitle(family_names[i]) +
        scale_color_manual(values = c("Endemic"="firebrick1",
                                   "Non endemic"="cadetblue"))+
        theme_bw()+
        theme(axis.title=element_blank(),
              axis.text=element_text(colour="black"),
              legend.title = element_text(size=8),
              legend.position = "inside",
              legend.position.inside = c(0.9,0.9),
              legend.box.background = element_blank())
    
    family_plots[[family_names[i]]] <- araneae_family
    ggsave(paste0("../figures/",family_names[i],"_distribution.png"), 
           plot=araneae_family, 
           height = 20, 
           width = 20,
           dpi = 300, 
           units="cm",
           device="png")

}


fig2 <- ggarrange(araneae_gr_base,
                  araneae_gr_gbif,
                  family_plots[["Gnaphosidae"]],
                  family_plots[["Dysderidae"]],
                  family_plots[["Theraphosidae"]],
                  family_plots[["Leptonetidae"]],
          labels = c("A", "B","C","D","E","F"),
          align = "hv",
          widths = c(1,1,1,1,1,1),
          ncol = 2,
          nrow = 3,
          font.label=list(color="black",size=22)) + bgcolor("white")

ggsave("../figures/Fig2.png", 
       plot=fig2, 
       height = 50, 
       width = 35,
       dpi = 300, 
       units="cm",
       device="png")

################################## region ###################################

## summary per region
region_area <- greece_regions |>
    mutate(area=round(set_units(st_area(geometry),km^2),0)) |>
    st_drop_geometry() |>
    distinct(NAME_2, area)

# sq_km
#
region_summary_1km <- araneae_gr_occ_tax |> 
    distinct(NAME_2,CELLCODE) |>
    group_by(NAME_2) |>
    summarise(sq_km=n())

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
    left_join(region_summary_end) |>
    left_join(region_summary_1km) |>
    left_join(region_area) |>
    mutate(sampled_proportion=round(sq_km/area,3)*100)

write_delim(region_all, "../results/regions_all_stats.tsv", delim="\t")

############################# References summary #########################

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

#################################### GBIF ####################################

gbif_cumulative <- gbif_araneae_gr |> 
    distinct(acceptedScientificName, collectionCode,year) |>
    arrange(year) |>
    mutate(Duplicates=duplicated(acceptedScientificName)) |> 
    mutate(First_occurrance=if_else(Duplicates=="FALSE",1,0)) |> 
    na.omit() |> 
    filter(First_occurrance==1) |>
    group_by(year) |> 
    summarise(species_year= n()) |>
    arrange(year) |> 
    mutate(Cumulative_occurrance= cumsum(species_year)) |> 
    mutate(Classification="GBIF species from specimen | sample | citation") |> 
    dplyr::select(-c(species_year)) |>
    distinct()

gbif_citizen_cumulative <- gbif_occ_citizen |> 
    distinct(acceptedScientificName, collectionCode,year) |>
    arrange(year) |>
    mutate(Duplicates=duplicated(acceptedScientificName)) |> 
    mutate(First_occurrance=if_else(Duplicates=="FALSE",1,0)) |> 
    na.omit() |> 
    filter(First_occurrance==1) |>
    group_by(year) |> 
    summarise(species_year= n()) |>
    arrange(year) |> 
    mutate(Cumulative_occurrance= cumsum(species_year)) |> 
    mutate(Classification="GBIF human observation species") |> 
    dplyr::select(-c(species_year)) |>
    distinct()

####################### araneae_accumulation ########################
araneae_accumulation <- bind_rows(araneae_ref_year,
                                  species_cumulative,
                                  endemic_cumulative_species,
                                  gbif_cumulative,
                                  gbif_citizen_cumulative)

############################ timeline figure ###############################
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
    scale_color_manual(values =c(
                                 "Endemic species to Greece"="firebrick1",
                                 "GBIF human observation species"="springgreen2",
                                 "GBIF species from specimen | sample | citation"="chartreuse4",
                                 "Publications"="darkorchid2",
                                 "All species"="gray0"
                                 ))+
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
          legend.position = c(0.23,0.87), 
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


