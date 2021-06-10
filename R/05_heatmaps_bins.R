#' ---
#' title: "Global Analysis of MPAs - Create climate summary heatmap"
#' author: "RS-eco"
#' ---

#Automatically install required packages, which are not yet installed
rm(list=ls()); invisible(gc())
packages <- c("sp", "raster", "tidyverse", "patchwork", "RStoolbox", "sf", "ggpmisc", "scico", "tagger")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load packages
l <- sapply(packages, require, character.only = TRUE, quietly=TRUE); rm(packages, l)

# Set working directory
workdir <- "/home/matt/Documents/futMPA"
setwd(workdir)

# Specify colour scheme
redblue <- scico(9, palette = 'roma')
redblue2 <- scico(255, palette = 'roma')
redwhiteblue <-  rev(scico(10, palette = 'vik'))
redwhiteblue2 <-  rev(scico(255, palette = 'vik'))

########################################

# Plot marine heatmap

# For the heatmap we need to summarise the data by temperature and salinity together
mar_dat <- readRDS("data/summary_all_marine_bins.rds")
colnames(mar_dat) <- c("temperature", "salinity", "I-II", "III-IV", "Not-designated", "Total", "V-VI", "n", "year", "rcp", "type")
mar_dat$marine <- "marine"

# Read and prepare data
mar_dat2 <- readRDS("data/summary_all_coast_bins.rds")
colnames(mar_dat2) <- c("temperature", "salinity", "I-II", "III-IV", "Not-designated", "Total", "V-VI", "n", "year", "rcp", "type")
mar_dat2$marine <- "coastal"

# Read and prepare data
mar_dat3 <- readRDS("data/summary_all_coastal+marine_bins.rds")
colnames(mar_dat3) <- c("temperature", "salinity", "I-II", "III-IV", "Not-designated", "Total", "V-VI", "n", "year", "rcp", "type")
mar_dat3$marine <- "coastal+marine"

# Merge data
mar_dat <- bind_rows(mar_dat, mar_dat2, mar_dat3) %>% unite("year_rcp", c(year, rcp), na.rm=T) %>%
    mutate(year_rcp = factor(year_rcp, levels=c("Present", "X2050AOGCM_RCP26", "X2050AOGCM_RCP60", "X2100AOGCM_RCP26", "X2100AOGCM_RCP60"), 
                             labels=c("2000-2014", "2040-2050 RCP2.6", "2040-2050 RCP6.0",
                                      "2090-2100 RCP2.6", "2090-2100 RCP6.0"))) #%>% drop_na()  

mar_dat$sum <- rowSums(mar_dat[,c("I-II", "III-IV", "V-VI", "Not-designated")], na.rm=T)

####################

# How can sum of individual protection be smaller than total?

mar_dat[round(mar_dat$sum, 1) < round(mar_dat$Total,1),]

####################

mar_dat$`I-II` <- ifelse(mar_dat$sum > mar_dat$Total, 
                         ifelse(mar_dat$`I-II` > mar_dat$Total, mar_dat$Total, mar_dat$`I-II`), 
                         mar_dat$`I-II`)
mar_dat$`III-IV` <- ifelse(mar_dat$sum > mar_dat$Total, 
                           ifelse(mar_dat[,c("I-II")] == mar_dat$Total, 0,
                                  ifelse(rowSums(mar_dat[,c("I-II", "III-IV")], na.rm=T) >= mar_dat$Total,
                                         mar_dat$Total-mar_dat$`I-II`,
                                         mar_dat$`III-IV`)), 
                           mar_dat$`III-IV`)
mar_dat$`V-VI` <- ifelse(mar_dat$sum > mar_dat$Total, 
                         ifelse(rowSums(mar_dat[,c("I-II", "III-IV")], na.rm=T) == mar_dat$Total, 0, 
                                ifelse(rowSums(mar_dat[,c("I-II", "III-IV", "V-VI")], na.rm=T) >= mar_dat$Total,
                                       mar_dat$Total-rowSums(mar_dat[,c("I-II", "III-IV")], na.rm=T),
                                       mar_dat$`V-VI`)), 
                         mar_dat$`V-VI`)
mar_dat$`Not-designated` <- ifelse(mar_dat$sum > mar_dat$Total, 
                                   ifelse(rowSums(mar_dat[,c("I-II", "III-IV", "V-VI")], na.rm=T) == mar_dat$Total, 0,
                                          ifelse(rowSums(mar_dat[,c("I-II", "III-IV", "V-VI", "Not-designated")], na.rm=T) >= mar_dat$Total,
                                                 mar_dat$Total-rowSums(mar_dat[,c("I-II", "III-IV", "V-VI")], na.rm=T),
                                                 mar_dat$`Not-designated`)), 
                                   mar_dat$`Not-designated`)

####################

# Add correct var values
df1 <- data.frame(salinity_x=seq(0, 41, by=1), salinity_y=seq(1, 42, by=1), salinity=1:42)
df2 <- data.frame(temperature_x=seq(-2, 33, by=1), temperature_y=seq(-1, 34, by=1), temperature=1:36)

mar_dat %<>% full_join(df1) %>% full_join(df2)
head(mar_dat)

####################

(no_cells <- mar_dat %>% group_by(year_rcp, type, marine) %>% summarise(no_cells=sum(n, na.rm=T)))

dat_sum <- mar_dat %>% group_by(salinity_x, salinity_y, temperature_x, temperature_y, year_rcp, type, marine) %>% 
    mutate(prot_cells = n*Total/100) %>% 
    summarise(cells=sum(n, na.rm=T), prot_cells = sum(prot_cells, na.rm=T)) %>%
    mutate(perc = prot_cells/cells*100) %>% mutate(perc = replace_na(perc, 0)) %>% 
    mutate(perc2 = as.character(cut(perc, breaks=c(0,1,2.5,5,10,20,30,50,100), include.lowest=T,
                                    labels=c("0 - 1", "1 - 2.5", "2.5 - 5", "5 - 10", "10 - 20",
                                             "20 - 30", "30 - 50", "50 - 100")))) %>% 
    mutate(perc2 = factor(if_else(perc == 0, "0", perc2),
                          levels=c("0", "0 - 1", "1 - 2.5", "2.5 - 5", "5 - 10", "10 - 20",
                                   "20 - 30", "30 - 50", "50 - 100"))) %>% drop_na()

p1 <- dat_sum %>% filter(year_rcp == "2000-2014", marine=="coastal+marine") %>% 
    ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, 
                             ymin=salinity_x, ymax=salinity_y, fill=perc2)) + 
    scale_fill_manual(name="% protected", values=redblue) + 
    facet_grid(type ~ .) + tag_facets(tag_pool = c("a", "c")) + 
    labs(x="Temperature", y="Salinity") + coord_cartesian() + theme_bw() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme(axis.title=element_text(size=12, face="bold"), 
          strip.text = element_text(size=12, face="bold"), strip.background = element_blank(), 
          panel.grid.minor = element_blank())

data <- dat_sum %>% filter(year_rcp == "2000-2014", marine=="coastal+marine", type=="Benthic") %>% 
    group_by(perc2) %>% summarise(area=sum(cells)*86005000)
data$fraction <- data$area / sum(data$area)
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2
data$label <- paste0(data$perc2, "\n", round(data$fraction*100,1))
p2 <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=perc2)) +
    geom_rect() + ggtitle("b) \t     Benthic area") + 
    # geom_text( x=2, aes(y=labelPosition, label=label, color=perc2), size=6) + # x here controls label position (inner / outer)
    scale_fill_manual(name="% protected", values=redblue) + coord_polar(theta="y") +
    xlim(c(1, 4)) + theme_void() + theme(legend.position = "none")
data <- dat_sum %>% filter(year_rcp == "2000-2014", marine=="coastal+marine", type=="Surface") %>% 
    group_by(type, marine, perc2) %>% summarise(area=sum(cells)*86005000)
data$fraction <- data$area / sum(data$area)
data$ymax <- cumsum(data$fraction)
data$ymin <- c(0, head(data$ymax, n=-1))
data$labelPosition <- (data$ymax + data$ymin) / 2
data$label <- paste0(data$perc2, "\n", round(data$fraction*100,1))
p3 <- ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=perc2)) +
    geom_rect() + ggtitle("d) \t     Surface area") + 
    # geom_text( x=2, aes(y=labelPosition, label=label, color=perc2), size=6) + # x here controls label position (inner / outer)
    scale_fill_manual(name="% protected", values=redblue) + coord_polar(theta="y") +
    xlim(c(1, 4)) + theme_void() + theme(legend.position = "none")

p1 + (p2 / p3)
ggsave(paste0("figures/Figure3_bins.png"), width=8, height=6, dpi=1000)

dat_sum %>% filter(marine=="coastal+marine") %>% 
    ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, 
                             ymin=salinity_x, ymax=salinity_y, fill=perc2)) + 
    scale_fill_manual(name="% protected", values=redblue) + 
    facet_grid(type ~ year_rcp) + tag_facets() + 
    labs(x="Temperature", y="Salinity") + coord_cartesian() + theme_bw() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme(axis.title=element_text(size=12, face="bold"), 
          strip.text = element_text(size=12, face="bold"), strip.background = element_blank(), 
          panel.grid.minor = element_blank())
ggsave(paste0("figures/fut_heatmap_coastal+marine_bins.png"), width=10, height=6, dpi=1000)

dat_sum %>% filter(marine=="coastal+marine") %>% 
    group_by(year_rcp, type, perc2) %>% summarise(area=sum(cells)*86005000) %>% drop_na() %>%
    ggplot(aes(x=year_rcp, y=(area/10000)/1e7, fill=perc2)) + 
    facet_grid(.~type, scales="free") + tag_facets() + 
    geom_bar(width=0.95, stat="identity", position=position_stack(reverse=T)) + 
    scale_fill_manual(name="% protected", values=redblue, 
                      guide = guide_legend(reverse = TRUE)) + 
    theme_bw() + labs(y="Area (10000 km2)") + 
    scale_x_discrete(expand=expansion(mult=c(.35,.35))) +
    scale_y_continuous(expand=expansion(mult=c(0,.01))) + theme_bw() + 
    theme(axis.text.x = element_text(angle=45, hjust=1),
          axis.text.y = element_text(angle=90, hjust=0.5),
          axis.title.x = element_blank(), strip.text = element_text(size=12, face="bold"), 
          strip.background = element_blank())
ggsave(paste0("figures/fut_area_bins.png"), width=8, height=6, dpi=1000)

change_sum <- dat_sum %>% ungroup() %>% filter(marine=="coastal+marine") %>% 
    select(-c(marine, prot_cells, prot_cells, cells, perc2)) %>% spread(year_rcp, perc) %>% 
    mutate_at(c("2000-2014", "2040-2050 RCP2.6", "2040-2050 RCP6.0", "2090-2100 RCP2.6", "2090-2100 RCP6.0"), replace_na, 0) %>%
    #mutate_at(c("2000-2014", "2040-2050 RCP2.6", "2040-2050 RCP6.0", "2090-2100 RCP2.6", "2090-2100 RCP6.0"), list(~ .*86005000)) %>%
    mutate_at(vars("2040-2050 RCP2.6", "2040-2050 RCP6.0", "2090-2100 RCP2.6", "2090-2100 RCP6.0"), list(~ . - `2000-2014`)) %>% 
    dplyr::select(-"2000-2014") %>% 
    pivot_longer(c("2040-2050 RCP2.6", "2040-2050 RCP6.0", "2090-2100 RCP2.6", "2090-2100 RCP6.0"), 
                 names_to="year_rcp", values_to="perc_change") # %>% mutate(perc_change = (area_change/3.675871e+14)/100)
summary(change_sum)    

p1 <- change_sum %>% filter(year_rcp %in% c("2090-2100 RCP2.6", "2090-2100 RCP6.0")) %>% ungroup() %>% 
    filter(type=="Benthic") %>% drop_na() %>% 
    ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, 
                             ymin=salinity_x, ymax=salinity_y, fill=perc_change)) + 
    scale_fill_gradientn(name="% change in \n% protected", colors=redwhiteblue2, limits=c(-100,100)) + 
    facet_grid(. ~ year_rcp) + tag_facets(tag_pool=c("a", "b")) + 
    labs(x="Temperature", y="Salinity") + coord_cartesian() + theme_bw() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme(axis.title.y=element_text(size=12, face="bold"), 
          axis.title.x= element_blank(), strip.text.x = element_text(size=12, face="bold"), 
          strip.background = element_blank(), panel.grid.minor = element_blank(),
          tagger.panel.tag.text = element_text(face="bold", size = 12))
p2 <- change_sum %>% filter(year_rcp %in% c("2090-2100 RCP2.6", "2090-2100 RCP6.0")) %>% ungroup() %>% 
    filter(type=="Surface") %>% drop_na() %>% 
    ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, 
                             ymin=salinity_x, ymax=salinity_y, fill=perc_change)) + 
    scale_fill_gradientn(name="% change in \n% protected", colors=redwhiteblue2, limits=c(-100,100)) + 
    facet_grid(. ~ year_rcp) + tag_facets(tag_pool=c("c", "d")) + 
    labs(x="Temperature", y="Salinity") + coord_cartesian() + theme_bw() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme(axis.title=element_text(size=12, face="bold"), 
          strip.text = element_blank(), strip.background = element_blank(), 
          panel.grid.minor = element_blank(),
          tagger.panel.tag.text = element_text(face="bold", size = 12))
(p1 / p2) + plot_layout(guides = 'collect')
ggsave(paste0("figures/Figure4_bins.png"), width=8, height=6, dpi=1000)

p1 <- change_sum %>% filter(year_rcp %in% c("2040-2050 RCP2.6", "2040-2050 RCP6.0")) %>% ungroup() %>% 
    filter(type=="Benthic") %>% drop_na() %>% 
    ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, 
                             ymin=salinity_x, ymax=salinity_y, fill=perc_change)) + 
    scale_fill_gradientn(name="% change in \n% protected", colors=redwhiteblue2, limits=c(-100,100)) + 
    facet_grid(. ~ year_rcp) + tag_facets(tag_pool=c("a", "b")) + 
    labs(x="Temperature", y="Salinity") + coord_cartesian() + theme_bw() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme(axis.title.y=element_text(size=12, face="bold"), 
          axis.title.x= element_blank(), strip.text.x = element_text(size=12, face="bold"), 
          strip.background = element_blank(), panel.grid.minor = element_blank(),
          tagger.panel.tag.text = element_text(face="bold", size = 12))
p2 <- change_sum %>% filter(year_rcp %in% c("2040-2050 RCP2.6", "2040-2050 RCP6.0")) %>% ungroup() %>% 
    filter(type=="Surface") %>% drop_na() %>% 
    ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, 
                             ymin=salinity_x, ymax=salinity_y, fill=perc_change)) + 
    scale_fill_gradientn(name="% change in \n% protected", colors=redwhiteblue2, limits=c(-100,100)) + 
    facet_grid(. ~ year_rcp) + tag_facets(tag_pool=c("c", "d")) + 
    labs(x="Temperature", y="Salinity") + coord_cartesian() + theme_bw() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme(axis.title=element_text(size=12, face="bold"), 
          strip.text = element_blank(), strip.background = element_blank(), 
          panel.grid.minor = element_blank(),
          tagger.panel.tag.text = element_text(face="bold", size = 12))
(p1 / p2) + plot_layout(guides = 'collect')
ggsave(paste0("figures/Figure4_bins_2040-2050.png"), width=8, height=6, dpi=1000)

########################################

# Plot marine heatmap of available climate space, protected climate space and proportion protected

dat_sum2 <- mar_dat %>% ungroup() %>% filter(type == "Surface", marine == "coastal+marine") %>%
    group_by(salinity_x, salinity_y, temperature_x, temperature_y, year_rcp) %>% 
    mutate(prot_cells = n*Total/100) %>% 
    summarise(cells=sum(n, na.rm=T), prot_cells = sum(prot_cells, na.rm=T)) %>%
    mutate(perc = prot_cells/cells*100) %>% mutate(perc = replace_na(perc, 0)) %>% 
    mutate(clim_space=cells*8600500/1e+6,
           prot_space=prot_cells*8600500/1e+6,
           prot_space2=perc*clim_space/100) %>% drop_na()

p1 <- dat_sum2 %>% ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, ymin=salinity_x, ymax=salinity_y, fill=clim_space/10000)) + 
    facet_grid(.~year_rcp) + scale_fill_gradientn(name="Area \n(10000 km2)", colors=redblue2, limits=c(0,NA)) + 
    labs(y="") + coord_cartesian() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme_bw() + theme(legend.key.height = unit(1.15, 'cm'), legend.title = element_text(vjust=0.85),
                       legend.background = element_blank(), strip.background = element_blank(),
                       strip.text=element_text(size=12, face="bold"), axis.title.x = element_blank())
p2 <- dat_sum2 %>% ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, ymin=salinity_x, ymax=salinity_y, fill=prot_space/10000)) + 
    facet_grid(.~year_rcp) + scale_fill_gradientn(name="Area \n(10000 km2)", colors=redblue2, limits=c(0,NA)) + 
    labs(y="Salinity") + coord_cartesian() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme_bw() + theme(legend.key.height = unit(1.15, 'cm'), legend.title = element_text(vjust=0.85),
                       legend.background = element_blank(), strip.text = element_blank(),
                       axis.title.y=element_text(size=12, face="bold"), axis.title.x = element_blank())
p3 <- dat_sum2 %>% ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, ymin=salinity_x, ymax=salinity_y, fill=perc)) + 
    facet_grid(.~year_rcp) + scale_fill_gradientn(name="% protected", colors=redblue2, limits=c(0,NA)) + 
    labs(x="Temperature", y="") + coord_cartesian() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme_bw() + theme(legend.key.height = unit(1.15, 'cm'), legend.title = element_text(vjust=0.85),
                       axis.title.x=element_text(size=12, face="bold"), 
                       legend.background = element_blank(), strip.text = element_blank(),
                       plot.title = element_text(hjust = 0.5, face="bold", size=12))

p1 / p2 / p3
ggsave(paste0("figures/heatmap_coastal+marine_surface_bins.png"), width=10, height=9, dpi=1000)

dat_sum3 <- mar_dat %>% ungroup() %>% filter(type == "Benthic", marine == "coastal+marine") %>%
    group_by(salinity_x, salinity_y, temperature_x, temperature_y, year_rcp, type, marine) %>% 
    mutate(prot_cells = n*Total/100) %>% 
    summarise(cells=sum(n, na.rm=T), prot_cells = sum(prot_cells, na.rm=T)) %>%
    mutate(perc = prot_cells/cells*100) %>% mutate(perc = replace_na(perc, 0)) %>% 
    mutate(clim_space=cells*8600500/1e+6,
           prot_space=prot_cells*8600500/1e+6,
           prot_space2=perc*clim_space/100) %>% drop_na()

p1 <- dat_sum3 %>% ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, ymin=salinity_x, ymax=salinity_y, fill=clim_space/10000)) + 
    facet_grid(.~year_rcp) + scale_fill_gradientn(name="Area \n(10000 km2)", colors=redblue2, limits=c(0,NA)) + 
    labs(y="") + coord_cartesian() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme_bw() + theme(legend.key.height = unit(1.15, 'cm'), legend.title = element_text(vjust=0.85),
                       legend.background = element_blank(), strip.background = element_blank(),
                       strip.text=element_text(size=12, face="bold"), axis.title.x = element_blank())
p2 <- dat_sum3 %>% ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, ymin=salinity_x, ymax=salinity_y, fill=prot_space/10000)) + 
    facet_grid(.~year_rcp) + scale_fill_gradientn(name="Area \n(10000 km2)", colors=redblue2, limits=c(0,NA)) + 
    labs(y="Salinity") + coord_cartesian() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme_bw() + theme(legend.key.height = unit(1.15, 'cm'), legend.title = element_text(vjust=0.85),
                       legend.background = element_blank(), strip.text = element_blank(),
                       axis.title.y=element_text(size=12, face="bold"), axis.title.x = element_blank())
p3 <- dat_sum3 %>% ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, ymin=salinity_x, ymax=salinity_y, fill=perc)) + 
    facet_grid(.~year_rcp) + scale_fill_gradientn(name="% protected", colors=redblue2, limits=c(0,NA)) + 
    labs(x="Temperature", y="") + coord_cartesian() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme_bw() + theme(legend.key.height = unit(1.15, 'cm'), legend.title = element_text(vjust=0.85),
                       axis.title.x=element_text(size=12, face="bold"), 
                       legend.background = element_blank(), strip.text = element_blank(),
                       plot.title = element_text(hjust = 0.5, face="bold", size=12))

p1 / p2 / p3
ggsave(paste0("figures/heatmap_coastal+marine_benthic_bins.png"), width=10, height=9, dpi=1000)

########################################

# Plot marine heatmap of  change in climate space, protected climate space and proportion protected

dat_sum4 <- mar_dat %>% 
    group_by(salinity_x, salinity_y, temperature_x, temperature_y, year_rcp, type, marine) %>% 
    mutate(prot_cells = n*Total/100) %>% 
    summarise(cells=sum(n, na.rm=T), prot_cells = sum(prot_cells, na.rm=T)) %>%
    mutate(perc = prot_cells/cells*100) %>% mutate(perc = replace_na(perc, 0)) %>% 
    mutate(clim_space=cells*8600500/1e+6,
           prot_space=prot_cells*8600500/1e+6,
           prot_space2=perc*clim_space/100) %>% drop_na()

clim_change <- dat_sum4 %>% ungroup() %>% filter(type=="Surface", marine=="coastal+marine") %>%
    dplyr::select(-c(type, marine, prot_cells, prot_space2, perc, cells, prot_space)) %>% spread(year_rcp, clim_space) %>% 
    mutate_at(vars("2040-2050 RCP2.6", "2040-2050 RCP6.0", "2090-2100 RCP2.6", "2090-2100 RCP6.0"), list(~ . - `2000-2014`)) %>% dplyr::select(-"2000-2014") %>% 
    pivot_longer(c("2040-2050 RCP2.6", "2040-2050 RCP6.0", "2090-2100 RCP2.6", "2090-2100 RCP6.0"), names_to="year_rcp", values_to="cells")
prot_change <- dat_sum4 %>% ungroup() %>% filter(type=="Surface", marine=="coastal+marine") %>%
    dplyr::select(-c(type, marine, prot_space2, perc, cells, prot_cells, clim_space)) %>% spread(year_rcp, prot_space) %>% 
    mutate_at(vars("2040-2050 RCP2.6", "2040-2050 RCP6.0", "2090-2100 RCP2.6", "2090-2100 RCP6.0"), list(~ . - `2000-2014`)) %>% dplyr::select(-"2000-2014") %>% 
    pivot_longer(c("2040-2050 RCP2.6", "2040-2050 RCP6.0", "2090-2100 RCP2.6", "2090-2100 RCP6.0"), names_to="year_rcp", values_to="cells")
perc_change <- dat_sum4 %>% ungroup() %>% filter(type=="Surface", marine=="coastal+marine") %>%
    dplyr::select(-c(type, marine, prot_space2, prot_space, cells, prot_cells, clim_space)) %>% spread(year_rcp, perc) %>%
    mutate_at(vars("2040-2050 RCP2.6", "2040-2050 RCP6.0", "2090-2100 RCP2.6", "2090-2100 RCP6.0"), list(~ . - `2000-2014`)) %>% dplyr::select(-"2000-2014") %>% 
    pivot_longer(c("2040-2050 RCP2.6", "2040-2050 RCP6.0", "2090-2100 RCP2.6", "2090-2100 RCP6.0"), names_to="year_rcp", values_to="cells")

p1 <- clim_change %>% ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, ymin=salinity_x, ymax=salinity_y, fill=cells/10000)) + 
    facet_grid(.~year_rcp) + scale_fill_gradientn(name="Area \n(10000 km2)", colors=redwhiteblue) + 
    labs(y="") + coord_cartesian() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme_bw() + theme(legend.key.height = unit(1.15, 'cm'), legend.title = element_text(vjust=0.85),
                       legend.background = element_blank(), strip.background = element_blank(),
                       strip.text=element_text(size=12, face="bold"), axis.title.x = element_blank())
p2 <- prot_change %>% ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, ymin=salinity_x, ymax=salinity_y, fill=cells/10000)) + 
    facet_grid(.~year_rcp) + scale_fill_gradientn(name="Area \n(10000 km2)", colors=redwhiteblue) + 
    labs(y="Salinity") + coord_cartesian() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme_bw() + theme(legend.key.height = unit(1.15, 'cm'), legend.title = element_text(vjust=0.85),
                       legend.background = element_blank(), strip.text = element_blank(),
                       axis.title.y=element_text(size=12, face="bold"), axis.title.x = element_blank())
p3 <- perc_change %>% ggplot() + geom_rect(aes(xmin=temperature_x, xmax=temperature_y, ymin=salinity_x, ymax=salinity_y, fill=cells)) + 
    facet_grid(.~year_rcp) + scale_fill_gradientn(name="% protected", colors=redwhiteblue) + 
    labs(x="Temperature", y="") + coord_cartesian() + 
    scale_x_continuous(expand=c(0,0)) + scale_y_continuous(expand=c(0,0)) + 
    theme_bw() + theme(legend.key.height = unit(1.15, 'cm'), legend.title = element_text(vjust=0.85),
                       axis.title.x=element_text(size=12, face="bold"), 
                       legend.background = element_blank(), strip.text = element_blank(),
                       plot.title = element_text(hjust = 0.5, face="bold", size=12))

p1 / p2 / p3
ggsave(paste0("figures/clim_change_cur_fut_prot.png"), width=10, height=9, dpi=1000)

########################################