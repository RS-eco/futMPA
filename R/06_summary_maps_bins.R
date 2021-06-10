#' ---
#' title: "Create environmental protection maps"
#' author: "RS-eco"
#' ---

rm(list=ls()); invisible(gc())

#Automatically install required packages, which are not yet installed
packages <- c("sp", "raster", "tidyr", "dplyr", "magrittr", "dtplyr",
              "ggplot2", "sf", "patchwork", "RStoolbox", "scico", "tagger")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load packages
l <- sapply(packages, require, character.only = TRUE, quietly=TRUE); rm(packages, l)

# Set working directory
workdir <- "/home/matt/Documents/futMPA"
setwd(workdir)

# Specify colour scheme
bluered <- rev(scico(255, palette = 'roma'))
redblue <- scico(9, palette = 'roma')

# Obtain world map
outline <- rgeos::gPolygonize(rgeos::gNode(as(rworldmap::getMap(resolution = "high"), "SpatialLines")))
outline <- rgeos::gUnaryUnion(outline)
outline <- sf::st_as_sf(outline)
outline_moll <- sf::st_transform(outline, "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

########################################

# Plot marine protection maps

# Read and prepare data
mar_dat_ind <- readRDS("data/summary_ind_marine_bins.rds")
head(mar_dat_ind)
tail(mar_dat_ind)
colnames(mar_dat_ind) <- c("year", "rcp", "type", "path", "var", "I-II", "III-IV", "Not-designated", "Total", "V-VI", "n")
mar_dat_ind$marine <- "marine"
unique(mar_dat_ind$path)
unique(mar_dat_ind$year)

# Read and prepare data
mar_dat_ind2 <- readRDS("data/summary_ind_coast_bins.rds")
head(mar_dat_ind2)
tail(mar_dat_ind2)
colnames(mar_dat_ind2) <- c("year", "rcp", "type", "path", "var", "I-II", "III-IV", "Not-designated", "Total", "V-VI", "n")
mar_dat_ind2$marine <- "coastal"
unique(mar_dat_ind2$path)
unique(mar_dat_ind2$year)

# Read and prepare data
mar_dat_ind3 <- readRDS("data/summary_ind_coastal+marine_bins.rds")
head(mar_dat_ind3)
tail(mar_dat_ind3)
colnames(mar_dat_ind3) <- c("year", "rcp", "type", "path", "var", "I-II", "III-IV", "Not-designated", "Total", "V-VI", "n")
mar_dat_ind3$marine <- "coastal+marine"
unique(mar_dat_ind3$path)
unique(mar_dat_ind3$year)

# Merge data
mar_dat_ind <- bind_rows(mar_dat_ind, mar_dat_ind2, mar_dat_ind3) %>% unite("year_rcp", c(year, rcp), na.rm=T)
rm(mar_dat_ind2, mar_dat_ind3)
mar_dat_ind$sum <- rowSums(mar_dat_ind[,c("I-II", "III-IV", "V-VI", "Not-designated")], na.rm=T)

####################

# How can sum of individual protection be smaller than total?

mar_dat_ind[round(mar_dat_ind$sum, 1) < round(mar_dat_ind$Total,1),]

####################

mar_dat_ind$`I-II` <- ifelse(mar_dat_ind$sum > mar_dat_ind$Total, 
                             ifelse(mar_dat_ind$`I-II` > mar_dat_ind$Total, mar_dat_ind$Total, mar_dat_ind$`I-II`), 
                             mar_dat_ind$`I-II`)
mar_dat_ind$`III-IV` <- ifelse(mar_dat_ind$sum > mar_dat_ind$Total, 
                               ifelse(mar_dat_ind[,c("I-II")] == mar_dat_ind$Total, 0,
                                      ifelse(rowSums(mar_dat_ind[,c("I-II", "III-IV")], na.rm=T) >= mar_dat_ind$Total,
                                             mar_dat_ind$Total-mar_dat_ind$`I-II`,
                                             mar_dat_ind$`III-IV`)), 
                               mar_dat_ind$`III-IV`)
mar_dat_ind$`V-VI` <- ifelse(mar_dat_ind$sum > mar_dat_ind$Total, 
                             ifelse(rowSums(mar_dat_ind[,c("I-II", "III-IV")], na.rm=T) == mar_dat_ind$Total, 0, 
                                    ifelse(rowSums(mar_dat_ind[,c("I-II", "III-IV", "V-VI")], na.rm=T) >= mar_dat_ind$Total,
                                           mar_dat_ind$Total-rowSums(mar_dat_ind[,c("I-II", "III-IV")], na.rm=T),
                                           mar_dat_ind$`V-VI`)), 
                             mar_dat_ind$`V-VI`)
mar_dat_ind$`Not-designated` <- ifelse(mar_dat_ind$sum > mar_dat_ind$Total, 
                                       ifelse(rowSums(mar_dat_ind[,c("I-II", "III-IV", "V-VI")], na.rm=T) == mar_dat_ind$Total, 0,
                                              ifelse(rowSums(mar_dat_ind[,c("I-II", "III-IV", "V-VI", "Not-designated")], na.rm=T) >= mar_dat_ind$Total,
                                                     mar_dat_ind$Total-rowSums(mar_dat_ind[,c("I-II", "III-IV", "V-VI")], na.rm=T),
                                                     mar_dat_ind$`Not-designated`)), 
                                       mar_dat_ind$`Not-designated`)

# Number of cells
mar_dat_ind %>% group_by(year_rcp, type, marine, path) %>% summarise(no_cells=sum(n))

# Area summary
mar_dat_ind %<>% dplyr::select(-c(Total, sum)) %>% 
  tidyr::gather(iucn_cat, perc, -c(year_rcp, type, marine, path, var, n)) %>% 
  mutate(iucn_cat = factor(iucn_cat, levels=c("Not-designated", "V-VI", "III-IV", "I-II"),
                           labels=c("Non-designated", "V-VI", "III-IV", "I-II"))) %>% drop_na()

mar_dat_ind %>% filter(path=="Current.Velocity", marine=="coastal", iucn_cat == "I-II") %>% 
  ungroup() %>% dplyr::select(perc) %>% summary()

# Summary
mar_dat_ind %>% group_by(path, year_rcp, type, marine, var) %>% summarise(sum=sum(perc)) %>% 
  summarise(max(sum))

# Plot marine heatmap

# For the heatmap we need to summarise the data by temperature and salinity together
mar_dat_all <- readRDS("data/summary_all_marine_bins.rds")
colnames(mar_dat_all) <- c("temperature", "salinity", "I-II", "III-IV", "Not-designated", "Total", "V-VI", "n", "year", "rcp", "type")
mar_dat_all$marine <- "marine"

# Read and prepare data
mar_dat_all2 <- readRDS("data/summary_all_coast_bins.rds")
colnames(mar_dat_all2) <- c("temperature", "salinity", "I-II", "III-IV", "Not-designated", "Total", "V-VI", "n", "year", "rcp", "type")
mar_dat_all2$marine <- "coastal"

# Read and prepare data
mar_dat_all3 <- readRDS("data/summary_all_coastal+marine_bins.rds")
colnames(mar_dat_all3) <- c("temperature", "salinity", "I-II", "III-IV", "Not-designated", "Total", "V-VI", "n", "year", "rcp", "type")
mar_dat_all3$marine <- "coastal+marine"

# Merge data
mar_dat_all <- bind_rows(mar_dat_all, mar_dat_all2, mar_dat_all3) %>% unite("year_rcp", c(year, rcp), na.rm=T)
rm(mar_dat_all2, mar_dat_all3)
mar_dat_all$sum <- rowSums(mar_dat_all[,c("I-II", "III-IV", "V-VI", "Not-designated")], na.rm=T)

####################

# How can sum of individual protection be smaller than total?

mar_dat_all[round(mar_dat_all$sum, 1) < round(mar_dat_all$Total,1),]

####################

mar_dat_all$`I-II` <- ifelse(mar_dat_all$sum > mar_dat_all$Total, 
                             ifelse(mar_dat_all$`I-II` > mar_dat_all$Total, mar_dat_all$Total, mar_dat_all$`I-II`), 
                             mar_dat_all$`I-II`)
mar_dat_all$`III-IV` <- ifelse(mar_dat_all$sum > mar_dat_all$Total, 
                               ifelse(mar_dat_all[,c("I-II")] == mar_dat_all$Total, 0,
                                      ifelse(rowSums(mar_dat_all[,c("I-II", "III-IV")], na.rm=T) >= mar_dat_all$Total,
                                             mar_dat_all$Total-mar_dat_all$`I-II`,
                                             mar_dat_all$`III-IV`)), 
                               mar_dat_all$`III-IV`)
mar_dat_all$`V-VI` <- ifelse(mar_dat_all$sum > mar_dat_all$Total, 
                             ifelse(rowSums(mar_dat_all[,c("I-II", "III-IV")], na.rm=T) == mar_dat_all$Total, 0, 
                                    ifelse(rowSums(mar_dat_all[,c("I-II", "III-IV", "V-VI")], na.rm=T) >= mar_dat_all$Total,
                                           mar_dat_all$Total-rowSums(mar_dat_all[,c("I-II", "III-IV")], na.rm=T),
                                           mar_dat_all$`V-VI`)), 
                             mar_dat_all$`V-VI`)
mar_dat_all$`Not-designated` <- ifelse(mar_dat_all$sum > mar_dat_all$Total, 
                                       ifelse(rowSums(mar_dat_all[,c("I-II", "III-IV", "V-VI")], na.rm=T) == mar_dat_all$Total, 0,
                                              ifelse(rowSums(mar_dat_all[,c("I-II", "III-IV", "V-VI", "Not-designated")], na.rm=T) >= mar_dat_all$Total,
                                                     mar_dat_all$Total-rowSums(mar_dat_all[,c("I-II", "III-IV", "V-VI")], na.rm=T),
                                                     mar_dat_all$`Not-designated`)), 
                                       mar_dat_all$`Not-designated`)

####################

# Add correct var values
df1 <- data.frame(salinity_x=seq(0, 41, by=1), salinity_y=seq(1, 42, by=1), salinity=1:42)
df2 <- data.frame(temperature_x=seq(-2, 33, by=1), temperature_y=seq(-1, 34, by=1), temperature=1:36)

mar_dat_all %<>% full_join(df1) %>% full_join(df2)
head(mar_dat_all)
summary(mar_dat_all)

#(no_cells <- mar_dat_all %>% ungroup() %>% summarise(no_cells=sum(n)) %>% as.numeric())

####################

vars <- c("Temperature", "Salinity", "Temperature.Salinity")
lapply(c("marine", "coastal", "coastal+marine"), function(marine){
  lapply(c("Benthic", "Surface"), function(type){
    lapply(c("Present", "2050AOGCM.RCP26", "2050AOGCM.RCP60", "2100AOGCM.RCP26", "2100AOGCM.RCP60"), function(year_rcp){
      files <- list.files("extdata", pattern=paste(year_rcp, type, sep="."), full.names=T)
      temp <- files[grepl(files, pattern="Temperature.Mean.*\\_bins.tif")]
      sal <- files[grepl(files, pattern="Salinity.Mean.*\\_bins.tif")]
      dat <- stack(temp, sal)
      if(year_rcp == "Present"){
        year.rcp <- year_rcp
      } else{
        year.rcp <- paste0("X", sub("[.]", "_", year_rcp))
      }
      lapply(1:length(vars), function(x){
        if(!file.exists(paste0("extdata/", paste(year_rcp, type, marine, vars[x], sep="."), "_protected_bins.tif"))){
          if(x == 1){
            dat_sub <- dat[[x]]
            dat_df <- mar_dat_ind %>% 
              filter(marine == marine, type == type, path=="Temperature",
                     year_rcp == year.rcp) %>% group_by(var) %>% mutate(prot_cells = n*perc/100) %>% 
              summarise(cells=sum(n, na.rm=T), prot_cells = sum(prot_cells, na.rm=T)) %>%
              mutate(Total = prot_cells/cells*100) %>% mutate(Total = replace_na(Total, 0))
          } else if(x == 2){
            dat_sub <- dat[[x]]
            dat_df <- mar_dat_ind %>% 
              filter(marine == marine, type == type, path=="Salinity",
                     year_rcp == year.rcp) %>% group_by(var) %>% mutate(prot_cells = n*perc/100) %>% 
              summarise(cells=sum(n, na.rm=T), prot_cells = sum(prot_cells, na.rm=T)) %>%
              mutate(Total = prot_cells/cells*100) %>% mutate(Total = replace_na(Total, 0))
          } else if(x == 3){
            dat_df <- mar_dat_all %>% filter(marine == marine, type == type,
                                             year_rcp == year.rcp) %>% 
              drop_na() %>% mutate(salinity = salinity*100) %>%
              mutate(var = temperature+salinity) %>% group_by(var) %>% 
              mutate(prot_cells = n*Total/100) %>% 
              summarise(cells=sum(n, na.rm=T), prot_cells = sum(prot_cells, na.rm=T)) %>%
              mutate(Total = prot_cells/cells*100) %>% mutate(Total = replace_na(Total, 0))
            dat_sub <- dat[[1]] + (dat[[2]]*100)
          }
          dat_df <- dat_df %>% ungroup() %>% dplyr::select(var, Total) %>% drop_na()
          summary(dat_df)
          
          # Subsitute values
          dat_sub <- rasterDT::subsDT(dat_sub, dat_df, by="var")
          
          # Save substituted rasters to file
          raster::writeRaster(dat_sub, 
                              filename=paste0("extdata/", paste(year_rcp, type, marine, vars[x], sep="."), "_protected_bins.tif"), 
                              format="GTiff", overwrite=T)
        }
      }); gc()
    })
  })
})

# Load substituted raster files
vars <- c("Temperature", "Salinity", "Temperature.Salinity")
prot_list <- lapply(c("marine", "coastal", "coastal+marine"), function(marine){
  dat <- lapply(c("Benthic", "Surface"), function(type){
    dat <- lapply(c("Present", "2050AOGCM.RCP26", "2050AOGCM.RCP60", "2100AOGCM.RCP26", "2100AOGCM.RCP60"), function(year_rcp){
      prot_all <- raster::stack(paste0("extdata/", paste(year_rcp, type, marine, vars, sep="."), "_protected_bins.tif")) %>% rasterToPoints() %>% as.data.frame(); gc()
      colnames(prot_all) <-c("x", "y", "Temperature", "Salinity", "Temperature + Salinity")
      prot_all$year_rcp <- year_rcp
      prot_all$type <- type
      prot_all$marine <- marine
      prot_all <- prot_all %>% tidyr::gather("var", "perc", -c(x,y,year_rcp, type, marine)) %>%
        mutate(var = factor(var, levels=c("Temperature", "Salinity", "Temperature + Salinity"))) %>%
        mutate(perc2 = as.character(cut(perc, breaks=c(0,1,2.5,5,10,20,30,50,100), include.lowest=T,
                                        labels=c("0 - 1", "1 - 2.5", "2.5 - 5", "5 - 10", "10 - 20",
                                                 "20 - 30", "30 - 50", "50 - 100")))) %>%
        mutate(perc2 = factor(if_else(perc == 0, "0", perc2),
                              levels=c("0", "0 - 1", "1 - 2.5", "2.5 - 5", "5 - 10", "10 - 20",
                                       "20 - 30", "30 - 50", "50 - 100")))
      return(prot_all)
    })
    bind_rows(dat)
  })
  bind_rows(dat)
}); gc()
prot_long <- bind_rows(prot_list); rm(prot_list); gc()
head(prot_long)
saveRDS(prot_long, "data/prot_maps_all.rds", compress="xz"); rm(prot_long); gc()

# Plot maps
prot_long <- readRDS("data/prot_maps_all.rds") %>% filter(marine == "coastal+marine", year_rcp=="2050AOGCM.RCP26")
p1 <- prot_long %>% ggplot() + geom_tile(aes(x=x, y=y, fill=perc2)) + 
  facet_grid(var~type, switch="y") + tag_facets() + 
  scale_fill_manual(values = redblue, name = "% Protected", na.value="transparent",
                    guide = guide_legend(reverse = TRUE)) + 
  geom_sf(data=outline_moll, fill=NA, color = "black", size = 1/.pt) +
  coord_sf() + theme_bw() + 
  theme(strip.background=element_blank(), axis.title=element_blank(), 
        strip.text = element_text(size=12, face="bold"), panel.border=element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank())
# Save figure
ggsave("figures/prot_maps_coastal+marine_2050-RCP26_bins.png", p1, width=9, height=6, dpi=600)

vars <- c("Temperature", "Salinity", "Temperature.Salinity")
prot_list <- lapply(c("marine", "coastal", "coastal+marine"), function(marine){
  dat <- lapply(c("Benthic", "Surface"), function(type){
    dat <- lapply(c("Present", "2050AOGCM.RCP26", "2050AOGCM.RCP60", "2100AOGCM.RCP26", "2100AOGCM.RCP60"), function(year_rcp){
      prot_all <- raster::stack(paste0("extdata/", paste(year_rcp, type, marine, vars, sep="."), "_protected_bins.tif")) %>% rasterToPoints() %>% as.data.frame(); gc()
      colnames(prot_all) <-c("x", "y", "Temperature", "Salinity", "Temperature + Salinity")
      prot_all$year_rcp <- year_rcp
      prot_all$type <- type
      prot_all$marine <- marine
      return(prot_all)
    })
    dat <- bind_rows(dat) 
    dat %>% tidyr::gather("var", "perc", -c(x,y,year_rcp, type, marine)) %>%
      mutate(year_rcp = factor(year_rcp, levels=c("Present", "2050AOGCM.RCP26", "2050AOGCM.RCP60", "2100AOGCM.RCP26", "2100AOGCM.RCP60"), 
                               labels=c("2000-2014", "2040-2050 RCP2.6", "2040-2050 RCP6.0",
                                        "2090-2100 RCP2.6", "2090-2100 RCP6.0"))) %>% 
      tidyr::pivot_wider(names_from=year_rcp, values_from=perc) %>%
      mutate_at(c("2000-2014", "2040-2050 RCP2.6", "2040-2050 RCP6.0", "2090-2100 RCP2.6", "2090-2100 RCP6.0"), replace_na, 0) %>%
      mutate_at(vars("2040-2050 RCP2.6", "2040-2050 RCP6.0", "2090-2100 RCP2.6", "2090-2100 RCP6.0"), list(~ . - `2000-2014`)) %>% 
      dplyr::select(-"2000-2014") %>% 
      pivot_longer(c("2040-2050 RCP2.6", "2040-2050 RCP6.0", "2090-2100 RCP2.6", "2090-2100 RCP6.0"), names_to="year_rcp", values_to="change")
    
  })
  bind_rows(dat)
}); gc()
prot_change <- bind_rows(prot_list); rm(prot_list); gc()
head(prot_change)
saveRDS(prot_change, "data/prot_change_maps_all.rds", compress="xz")

prot_change <- readRDS("data/prot_change_maps_all.rds")
prot_change %<>% 
  mutate(var = factor(var, levels=c("Temperature", "Salinity", "Temperature + Salinity"))) %>%
  mutate(perc2 = factor(cut(change, breaks=c(-60, -20, -10, -5, 0,5,10,20,61), include.lowest=T,
                                  labels=c("-60 - -20", "-20 - -10", "-10 - -5", "-5 - 0", "0 - 5", "5 - 10",
                                           "10 - 20", "20 - 60")), 
                          levels=c("-60 - -20", "-20 - -10", "-10 - -5", "-5 - 0", "0 - 5", "5 - 10",
                                   "10 - 20", "20 - 60")))

p1 <- prot_change %>% filter(year_rcp == "2040-2050 RCP2.6") %>% 
  ggplot() + geom_tile(aes(x=x, y=y, fill=perc2)) + 
  facet_grid(var~type, switch="y") + tag_facets() + 
  scale_fill_manual(values = redblue, name = "% Protected", na.value="transparent",
                    guide = guide_legend(reverse = TRUE)) + 
  geom_sf(data=outline_moll, fill=NA, color = "black", size = 1/.pt) +
  coord_sf() + theme_bw() + 
  theme(strip.background=element_blank(), axis.title=element_blank(), 
        strip.text = element_text(size=12, face="bold"), panel.border=element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank())

# Save figure
ggsave("figures/Figure5_bins.png", p1, width=9, height=6, dpi=600)
