#' ---
#' title: "Global Analysis of MPAs - Create environmental summary figures"
#' author: "RS-eco"
#' ---

rm(list=ls()); invisible(gc())

#Automatically install required packages, which are not yet installed
#devtools::install_github("eliocamp/tagger")
packages <- c("tidyverse", "patchwork", "ggpubr", "ggpmisc", "magrittr", "tagger")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load packages
l <- sapply(packages, require, character.only = TRUE, quietly=TRUE); rm(packages, l)

# Set working directory
workdir <- "/home/matt/Documents/futMPA"
setwd(workdir)

col5 <- c('black', "#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF")
col4 <- c("#00468BFF", "#ED0000FF", "#42B540FF", "#0099B4FF")
col3 <- c("black", "#ED0000FF", "#0099B4FF")
col2 <- c("black", "red")

########################################

# Plot of marine data

# Read and prepare data
mar_dat <- readRDS("data/summary_ind_marine_bins.rds")
head(mar_dat)
tail(mar_dat)
colnames(mar_dat) <- c("year", "rcp", "type", "path", "var", "I-II", "III-IV", "Not-designated", "Total", "V-VI", "n")
mar_dat$marine <- "marine"
unique(mar_dat$path)
unique(mar_dat$year)

# Read and prepare data
mar_dat2 <- readRDS("data/summary_ind_coast_bins.rds")
head(mar_dat2)
tail(mar_dat2)
colnames(mar_dat2) <- c("year", "rcp", "type", "path", "var", "I-II", "III-IV", "Not-designated", "Total", "V-VI", "n")
mar_dat2$marine <- "coastal"
unique(mar_dat2$path)
unique(mar_dat2$year)

# Read and prepare data
mar_dat3 <- readRDS("data/summary_ind_coastal+marine_bins.rds")
head(mar_dat3)
tail(mar_dat3)
colnames(mar_dat3) <- c("year", "rcp", "type", "path", "var", "I-II", "III-IV", "Not-designated", "Total", "V-VI", "n")
mar_dat3$marine <- "coastal+marine"
unique(mar_dat3$path)
unique(mar_dat3$year)

# Merge data
mar_dat <- bind_rows(mar_dat, mar_dat2, mar_dat3) %>% unite("year_rcp", c(year, rcp), na.rm=T) %>%
  mutate(year_rcp = factor(year_rcp, levels=c("Present", "X2050AOGCM_RCP26", "X2050AOGCM_RCP60", "X2100AOGCM_RCP26", "X2100AOGCM_RCP60"), 
                           labels=c("2000-2014", "2040-2050 RCP2.6", "2040-2050 RCP6.0", "2090-2100 RCP2.6", "2090-2100 RCP6.0"))) #%>% drop_na()  
rm(mar_dat2, mar_dat3); gc()
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

# Number of cells
mar_dat %>% group_by(year_rcp, type, marine, path) %>% summarise(no_cells=sum(n))

# Area summary
mar_dat %<>% dplyr::select(-c(Total, sum)) %>% 
  tidyr::gather(iucn_cat, perc, -c(year_rcp, type, marine, path, var, n)) %>% 
  mutate(iucn_cat = factor(iucn_cat, levels=rev(c("I-II", "III-IV",  "V-VI", "Not-designated")),
                           labels=rev(c("I-II", "III-IV",  "V-VI", "Not-designated")))) %>% drop_na()

mar_dat %>% filter(path=="Current.Velocity", marine=="coastal+marine", iucn_cat == "I-II") %>% 
  ungroup() %>% dplyr::select(perc) %>% summary()

# Summary
mar_dat %>% group_by(path, year_rcp, type, marine, var) %>% summarise(sum=sum(perc)) %>% 
  summarise(max(sum))

####################

# Calculate % area coverage per bin

# Total number of cells
(tot_sum <- mar_dat %>% group_by(path, year_rcp, type, marine, var) %>% summarise(n=sum(n)/4) %>%
   group_by(path, year_rcp, type, marine) %>% summarise(sum=sum(n)))
max(tot_sum$sum)
min(tot_sum$sum)

# % number of cells per bin
(area_val <- mar_dat %>% group_by(path, year_rcp, type, marine, var) %>% summarise(n=sum(n)/4) %>% 
    left_join(tot_sum) %>% mutate(perc_area=(n/sum)*100))

# Test if area adds up to 100 %
area_val %>% summarise(sum(perc_area))

# Area per cell
(area_cell <- 86005000) #m2

# Number of cells
(number_cells <- max(tot_sum$sum))

# Total area (m2)
(tot_area <- area_cell*number_cells)

# Total area (km2)
tot_area/10000000

# Total area (10000000 km2)
tot_area1e7 <- tot_area/(1000000*10000000)

# Total area in 10000 km2
area_val %>% summarise(sum(perc_area*tot_area1e7))

####################

#' ### Hyper-geometric distribution

# Total area & Total area protected
(tot_sum <- mar_dat %>% group_by(path, year_rcp, type, marine, var) %>% summarise(sum=sum(n)/n(), prot_cells=sum((perc*n)/100)) %>%
   ungroup() %>% group_by(path, year_rcp, type, marine) %>% summarise(sum=sum(sum), prot_cells=sum(prot_cells)) %>%
   mutate(prop_prot=(prot_cells/sum)*100))

# Global area climate
(clim_area <- mar_dat %>% group_by(path, year_rcp, type, marine, var) %>% summarise(area_clim=sum(n)/n()) %>% 
    left_join(tot_sum))

# Expected proportion climate  
(exp_val <- clim_area %>% mutate(perc_clim = area_clim/sum*100) %>% 
    mutate(exp = (perc_clim*prop_prot),
           exp_aichi = 0.15*perc_clim,
           var_exp = (perc_clim*prop_prot*(1-perc_clim)*(1-prop_prot)/(sum-1))))

# Need to multiply by 100 to get perc value rather than proportion

# => Need to consider number of bins, not just 100!!!

# Need to multiply by bin size to get values up to 100 %
exp_val %>% group_by(path) %>% summarise(sum(perc_clim))
exp_val %>% group_by(path) %>% summarise(sum(exp))
exp_val %>% group_by(path) %>% summarise(sum(exp_aichi))

####################

# Add correct var values
vars <- c("Chlorophyll", "Cloud.cover", "Current.Velocity", "Diffuse.attenuation", "Dissolved.oxygen", "Ice.cover", "Ice.thickness", "Iron",                
          "Light.bottom", "Nitrate", "Phosphate", "Phytoplankton", "Primary.productivity", "Salinity", "Silicate", "Temperature")       
m <- list(data.frame(x=seq(0, 6.75, by=0.25), y=seq(0.25, 7, by=0.25), z=1:28),
          data.frame(x=seq(-5, 1.75, by=0.25), y=seq(-4.75, 2, by=0.25), z=1:28),
          data.frame(x=seq(0, 1.9, by=0.1), y=seq(0.1, 2, by=0.1), z=1:20),
          data.frame(x=seq(-1, 17.5, by=0.5), y=seq(-0.5, 18, by=0.5), z=1:38),
          data.frame(x=seq(0, 413, by=7), y=seq(7, 420, by=7), z=1:60),
          data.frame(x=seq(-1, 0.9, by=0.1), y=seq(-0.9, 1, by=0.1), z=1:20),
          data.frame(x=seq(-1, 7.5, by=0.5), y=seq(-0.5, 8, by=0.5), z=1:18),
          data.frame(x=seq(0, 0.95, by=0.05), y=seq(0.05, 1, by=0.05), z=1:20),
          data.frame(x=seq(-1, 53, by=1), y=seq(0, 54, by=1), z=1:55),
          data.frame(x=seq(-1, 107, by=2), y=seq(1, 109, by=2), z=1:55),
          data.frame(x=seq(0, 3.8, by=0.2), y=seq(0.2, 4, by=0.2), z=1:20),
          data.frame(x=seq(0, 19, by=1), y=seq(1, 20, by=1), z=1:20),
          data.frame(x=seq(-1, 0.9, by=0.1), y=seq(-0.9, 1, by=0.1), z=1:20),
          data.frame(x=seq(0, 41, by=1), y=seq(1, 42, by=1), z=1:42),
          data.frame(x=seq(0, 310, by=10), y=seq(10, 320, by=10), z=1:32),
          data.frame(x=seq(-2, 33, by=1), y=seq(-1, 34, by=1), z=1:36))
names(m) <- vars
df <- bind_rows(m, .id="path")
colnames(df) <- c("path", "x", "y", "var")

mar_dat %<>% left_join(df) %>% mutate(var2=rowMeans(cbind(x,y))) %>% full_join(area_val)
head(mar_dat)
levels(mar_dat$iucn_cat)

####################

# Plot summary figures

sub_dat <- mar_dat %>% filter(year_rcp=="2000-2014", path == "Temperature", type=="Surface", marine == "coastal+marine")
sub_dat_group <- mar_dat %>% filter(year_rcp=="2000-2014", path == "Temperature", type=="Surface", marine == "coastal+marine") %>%
  group_by(var2) %>% summarise(perc_area=sum(perc_area)/n(), perc=sum(perc))
#setDT(sub_dat_group)[, delta_perc := perc - lag(perc, 1L)]
#setDT(sub_dat_group)[, delta_var2 := var2 - lag(var2, 1L)]
p1 <- ggplot() + geom_bar(data=sub_dat_group, aes(x=var2, y=perc_area*tot_area1e7), 
                    colour="black", fill="lightgrey", stat="identity", width=1) +
  geom_bar(data=sub_dat, aes(x = var2, y=perc*(perc_area/100)*tot_area1e7, fill=iucn_cat), 
           stat="identity", position = "stack", width=1) + 
  geom_bar(data=sub_dat_group, aes(x=var2, y=perc*(perc_area/100)*tot_area1e7), 
           col="black", fill="transparent", stat="identity", width=1) +
  geom_point(data=sub_dat_group, aes(x=var2, y=perc*tot_area1e7), col="red") +
  geom_line(data=sub_dat_group, aes(x=var2, y=perc*tot_area1e7), col="red", lty="dashed") + 
  #geom_segment(data=sub_dat_group, aes(x=(var2 - delta_var2), y=(perc - delta_perc)*tot_area1e7,
  #                                     xend=var2, yend=(perc)*tot_area1e7), col="red") +
  xlab("") + scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), 
                     sec.axis = sec_axis(~ ./tot_area1e7, breaks=c(0,2.5,5,7.5,10,12.5,15, 17.5), name="% protected")) + 
  scale_fill_manual(name="IUCN", values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"),
                    guide=guide_legend(reverse=T)) + 
  theme_bw() + theme(legend.position="bottom", axis.title.y=element_blank(), 
                     axis.title.x=element_text(size=12, face="bold"), 
                     panel.grid.minor = element_blank())
p1

# Create separate Figure legend
leg <- ggpubr::as_ggplot(ggpubr::get_legend(p1))

# Create figure for each variable

####

# Iron and Light.bottom currently only have one value!!!
# Primary.productivity only has 2 values!!!
# Cloud.cover not informative

# => Have all been removed

###

vars <- c("Chlorophyll", "Diffuse.attenuation", "Dissolved.oxygen", "Ice.cover", "Current.Velocity", "Ice.thickness",
          "Nitrate", "Phosphate", "Phytoplankton", "Silicate", "Temperature", "Salinity")
p <- lapply(1:length(vars), function(i){
  let <- letters[i]
  sub_dat <- mar_dat %>% filter(year_rcp=="2000-2014", marine == "coastal+marine", type=="Surface", path == vars[i]) 
  wid <- sub_dat$var2[2] - sub_dat$var2[1]
  sub_dat_group <- sub_dat %>% group_by(var2) %>% summarise(perc_area=sum(perc_area)/n(), perc=sum(perc))
  p <- ggplot() + geom_bar(data=sub_dat_group, aes(x=var2, y=perc_area*tot_area1e7), 
                            colour="black", fill="lightgrey", stat="identity", width=wid) +
    geom_bar(data=sub_dat, aes(x = var2, y=perc*(perc_area/100)*tot_area1e7, fill=iucn_cat), 
             stat="identity", position = "stack", width=wid) + 
    geom_bar(data=sub_dat_group, aes(x=var2, y=perc*(perc_area/100)*tot_area1e7), 
             col="black", fill="transparent", stat="identity", width=wid) +
    geom_point(data=sub_dat_group, aes(x=var2, y=perc*tot_area1e7), col="red") +
    geom_line(data=sub_dat_group, aes(x=var2, y=perc*tot_area1e7), col="red", lty="dashed") + 
    geom_text_npc(aes(npcx = 0.05, npcy=0.95, label=paste0(let, ")"))) + 
    xlab(vars[i]) + scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=expansion(mult=c(0,.01)), 
                       sec.axis = sec_axis(~ ./tot_area1e7)) + 
    scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"),
                      guide=guide_legend(reverse=T)) + 
    theme_bw() + theme(legend.position = "none", axis.title.y=element_blank(), 
                       axis.title.x=element_text(size=12, face="bold"), 
                       panel.grid.minor = element_blank()) + coord_trans(y="sqrt")
  if(i == 5){
    p <- p + ylab("Area (1e7 km2)") + theme(axis.title.y=element_text(face="bold", size=12, angle = 90))
  } else if(i == 8){
    p <- p + scale_y_continuous(expand=expansion(mult=c(0,.01)), 
                                sec.axis = sec_axis(~ ./tot_area1e7, name="Protected (%)")) + 
      ylab("") + theme(axis.title.y=element_text(face="bold", size=12, angle = 90))
  }
  return(p)
})
p <- (wrap_plots(p) / leg) + plot_layout(heights=c(9,1)) #+ plot_annotation(tag_levels = 'a')
p
# source("R/add_global_label.R")
# p %>% add_global_label(Ylab ="\t Area (1e07 km2)", 
ggsave("figures/Figure1_bins.png", dpi=1000, width=10, height=6)

vars <- c("Chlorophyll", "Dissolved.oxygen", "Current.Velocity", 
          "Nitrate", "Phosphate", "Phytoplankton", "Silicate", "Temperature", "Salinity")
p <- lapply(1:length(vars), function(i){
  let <- letters[i]
  sub_dat <- mar_dat %>% filter(year_rcp=="2000-2014", marine == "coastal+marine", type=="Benthic", path == vars[i]) 
  wid <- sub_dat$var2[2] - sub_dat$var2[1]
  sub_dat_group <- sub_dat %>% group_by(var2) %>% summarise(perc_area=sum(perc_area)/n(), perc=sum(perc))
  p <- ggplot() + geom_bar(data=sub_dat_group, aes(x=var2, y=perc_area*tot_area1e7), 
                           colour="black", fill="lightgrey", stat="identity", width=wid) +
    geom_bar(data=sub_dat, aes(x = var2, y=perc*(perc_area/100)*tot_area1e7, fill=iucn_cat), 
             stat="identity", position = "stack", width=wid) + 
    geom_bar(data=sub_dat_group, aes(x=var2, y=perc*(perc_area/100)*tot_area1e7), 
             col="black", fill="transparent", stat="identity", width=wid) +
    geom_point(data=sub_dat_group, aes(x=var2, y=perc*tot_area1e7), col="red") +
    geom_line(data=sub_dat_group, aes(x=var2, y=perc*tot_area1e7), col="red", lty="dashed") + 
    geom_text_npc(aes(npcx = 0.05, npcy=0.95, label=paste0(let, ")"))) + 
    xlab(vars[i]) + scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=expansion(mult=c(0,.01)), 
                       sec.axis = sec_axis(~ ./tot_area1e7)) + 
    scale_fill_manual(values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"),
                      guide=guide_legend(reverse=T)) + 
    theme_bw() + theme(legend.position = "none", axis.title.y=element_blank(), 
                       axis.title.x=element_text(size=12, face="bold"), 
                       panel.grid.minor = element_blank()) + coord_trans(y="sqrt")
  if(i == 4){
    p <- p + ylab("Area (1e7 km2)") + theme(axis.title.y=element_text(face="bold", size=12, angle = 90))
  } else if(i == 6){
    p <- p + scale_y_continuous(expand=expansion(mult=c(0,.01)), 
                                sec.axis = sec_axis(~ ./tot_area1e7, name="Protected (%)")) + 
      ylab("") + theme(axis.title.y=element_text(face="bold", size=12, angle = 90))
  }
  return(p)
})
p <- (wrap_plots(p) / leg) + plot_layout(heights=c(9,1)) #+ plot_annotation(tag_levels = 'a')
p
# source("R/add_global_label.R")
# p %>% add_global_label(Ylab ="\t Area (1e07 km2)", 
ggsave("figures/Figure1_bins_benthic.png", dpi=1000, width=10, height=6)

mar_dat$type <- factor(mar_dat$type)
str(mar_dat)
p1 <- mar_dat %>% filter(path == "Temperature", marine == "coastal+marine") %>% 
  separate(year_rcp, into=c("year", "rcp"), sep=" ") %>% 
  mutate(rcp = replace_na(rcp, "RCP2.6")) %>%
  group_by(year, rcp, path, marine, type, var2) %>% summarise(perc=sum(perc)) %>% 
  ggplot(aes(x = var2, y=perc, colour = year, lty=rcp)) + geom_line() + 
  facet_grid(path ~ type, scales="free") + 
  scale_y_continuous(expand=expansion(mult=c(0,.025))) + 
  scale_colour_manual(values=col3) + 
  theme_bw() + theme(legend.position = "none", strip.text = element_text(size=12, face="bold"), 
                     axis.title = element_blank(),
                     strip.background = element_blank(), strip.placement = "outside", 
                     panel.grid.minor = element_blank())
p2 <- mar_dat %>% filter(path == "Salinity", marine == "coastal+marine") %>% 
  separate(year_rcp, into=c("year", "rcp"), sep=" ") %>% 
  mutate(rcp = replace_na(rcp, "RCP2.6")) %>%
  group_by(year, rcp, path, marine, type, var2) %>% summarise(perc=sum(perc)) %>% 
  ggplot(aes(x = var2, y=perc, colour = year, lty=rcp)) + geom_line() + 
  facet_grid(path ~ type, scales="free") + 
  scale_y_continuous(expand=expansion(mult=c(0,.025))) + 
  scale_colour_manual(values=col3) + 
  theme_bw() + theme(legend.position = "none", strip.text.y = element_text(size=12, face="bold"), 
                     strip.text.x = element_blank(), 
                     axis.title = element_blank(), strip.background = element_blank(),
                     strip.placement = "outside", panel.grid.minor = element_blank())
p3 <- mar_dat %>% filter(path == "Current.Velocity", marine == "coastal+marine") %>% 
  separate(year_rcp, into=c("year", "rcp"), sep=" ") %>% 
  mutate(rcp = replace_na(rcp, "RCP2.6")) %>%
  group_by(year, rcp, path, marine, type, var2) %>% summarise(perc=sum(perc)) %>% 
  ggplot(aes(x = var2, y=perc, colour = year, lty=rcp)) + geom_line() + 
  facet_grid(path ~ type, scales="free") + 
  scale_y_continuous(expand=expansion(mult=c(0,.025))) + 
  scale_colour_manual(values=col3) + 
  theme_bw() + theme(legend.position = "none", strip.text.y = element_text(size=12, face="bold"), 
                     strip.text.x = element_blank(), 
                     axis.title = element_blank(), strip.background = element_blank(),
                     strip.placement = "outside", panel.grid.minor = element_blank())
p4 <- mar_dat %>% filter(path == "Ice.thickness", marine == "coastal+marine") %>% 
  separate(year_rcp, into=c("year", "rcp"), sep=" ") %>% 
  mutate(rcp = replace_na(rcp, "RCP2.6")) %>%
  group_by(year, rcp, path, marine, type, var2) %>% summarise(perc=sum(perc)) %>% 
  ggplot(aes(x = var2, y=perc, colour = year, lty=rcp)) + geom_line() + 
  facet_grid(path ~ type, scales="free", drop=F) + 
  scale_y_continuous(expand=expansion(mult=c(0,.025))) + 
  scale_colour_manual(name="Year", values=col3) + 
  labs(lty="RCP") + 
  theme_bw() + theme(legend.position = "bottom", strip.text.y = element_text(size=12, face="bold"), 
                     strip.text.x = element_blank(), 
                     axis.title = element_blank(), strip.background = element_blank(),
                     strip.placement = "outside", panel.grid.minor = element_blank())
leg2 <- ggpubr::get_legend(p4)
p <- (p1 / p2 / p3 / {p4 + theme(legend.position="none")}) / leg2 + plot_layout(heights=c(5,5,5,5,1))
source("R/add_global_label.R")
p %>% add_global_label(Ylab ="% area protected", size=12, Ygap=0.04) 
ggsave("figures/pres_fut_line_bins.png", dpi=1000, width=6, height=9)

mar_dat$type <- factor(mar_dat$type)
mar_dat2 <- mar_dat %>% ungroup() %>% select(-perc) %>%
  spread(year_rcp, n) %>% mutate_at(c("2000-2014", "2040-2050 RCP2.6", "2040-2050 RCP6.0",  "2090-2100 RCP2.6", "2090-2100 RCP6.0"), replace_na, 0) %>%
  mutate_at(vars("2040-2050 RCP2.6", "2040-2050 RCP6.0",  "2090-2100 RCP2.6", "2090-2100 RCP6.0"), funs(. - `2000-2014`)) %>% select(-"2000-2014") %>% 
  pivot_longer(c("2040-2050 RCP2.6", "2040-2050 RCP6.0",  "2090-2100 RCP2.6", "2090-2100 RCP6.0"), names_to="year_rcp", values_to="n")

p1 <- mar_dat2 %>% filter(path == "Current.Velocity", marine == "coastal+marine") %>% 
  separate(year_rcp, into=c("year", "rcp"), sep=" ") %>% 
  group_by(path, year, rcp, type, var2) %>% summarise(perc=sum(n)*86005000/3.675871e+14*100) %>% 
  ggplot(aes(x = var2, y=perc, colour = year, lty=rcp)) + geom_line() + 
  geom_hline(yintercept=0, colour="black") + 
  facet_grid(path ~ type, scales="free") + 
  tag_facets() + 
  scale_y_continuous(expand=expansion(mult=c(0,.025))) + 
  #, sec.axis = sec_axis(~ .*36758.71, name="Area (1e07 km2)")) + 
  scale_colour_manual(name="Year", values=col4[c(2,4)]) + 
  labs(lty="RCP") + 
  theme_bw() + theme(legend.position = "none", strip.text = element_text(size=12, face="bold"), 
                     axis.title = element_blank(), 
                     strip.background = element_blank(), strip.placement = "outside", 
                     panel.grid.minor = element_blank())
p2 <- mar_dat2 %>% filter(path == "Ice.thickness", marine == "coastal+marine") %>% 
  separate(year_rcp, into=c("year", "rcp"), sep=" ") %>% 
  group_by(path, year, rcp, type, var2) %>% summarise(perc=sum(n)*86005000/3.675871e+14*100) %>% 
  ggplot(aes(x = var2, y=perc, colour = year, lty=rcp)) + geom_line() + 
  geom_hline(yintercept=0, colour="black") + 
  facet_grid(path ~ type, scales="free", drop=F) + 
  tag_facets(tag_pool=c("c", "d")) + 
  scale_y_continuous(expand=expansion(mult=c(0,.025))) + 
  #, sec.axis = sec_axis(~ .*36758.71, name="Area (1e07 km2)")) + 
  scale_colour_manual(name="Year", values=col4[c(2,4)]) + 
  labs(lty="RCP") + 
  theme_bw() + theme(legend.position = "bottom", strip.text.y = element_text(size=12, face="bold"), 
                     axis.title = element_blank(), strip.text.x = element_blank(),
                     strip.background = element_blank(), strip.placement = "outside", 
                     panel.grid.minor = element_blank(),
                     tagger.panel.tag.text = element_text(face="bold", size = 12))
p <- p1 / p2
source("R/add_global_label.R")
p %>% add_global_label(Ylab ="% change in protected area", size=10, Ygap=0.04) 
ggsave("figures/perc_change_line_bins.png", dpi=1000, width=8, height=8)

p1 <- mar_dat2 %>% filter(path == "Temperature", marine == "coastal+marine") %>% 
  separate(year_rcp, into=c("year", "rcp"), sep=" ") %>% 
  group_by(path, year, rcp, type, var2) %>% summarise(perc=sum(n)*86005000/3.675871e+14*100) %>% 
  ggplot(aes(x = var2, y=perc, colour = year, lty=rcp)) + geom_line() + 
  geom_hline(yintercept=0, colour="black") + 
  facet_grid(path ~ type, scales="free") + 
  tag_facets() + 
  scale_y_continuous(expand=expansion(mult=c(0,.025))) + 
  #, sec.axis = sec_axis(~ .*36758.71, name="Area (1e07 km2)")) + 
  scale_colour_manual(name="Year", values=col4[c(2,4)]) + 
  labs(lty="RCP") + 
  theme_bw() + theme(legend.position = "none", strip.text = element_text(size=12, face="bold"), 
                     axis.title = element_blank(), 
                     strip.background = element_blank(), strip.placement = "outside", 
                     panel.grid.minor = element_blank())
p2 <- mar_dat2 %>% filter(path == "Salinity", marine == "coastal+marine") %>% 
  separate(year_rcp, into=c("year", "rcp"), sep=" ") %>% 
  group_by(path, year, rcp, type, var2) %>% summarise(perc=sum(n)*86005000/3.675871e+14*100) %>% 
  ggplot(aes(x = var2, y=perc, colour = year, lty=rcp)) + geom_line() + 
  geom_hline(yintercept=0, colour="black") + 
  facet_grid(path ~ type, scales="free", drop=F) + 
  tag_facets(tag_pool=c("c", "d")) + 
  scale_y_continuous(expand=expansion(mult=c(0,.025))) + 
  #, sec.axis = sec_axis(~ .*36758.71, name="Area (1e07 km2)")) + 
  scale_colour_manual(name="Year", values=col4[c(2,4)]) + 
  labs(lty="RCP") + 
  theme_bw() + theme(legend.position = "bottom", strip.text.y = element_text(size=12, face="bold"), 
                     axis.title = element_blank(), strip.text.x = element_blank(),
                     strip.background = element_blank(), strip.placement = "outside", 
                     panel.grid.minor = element_blank(),
                     tagger.panel.tag.text = element_text(face="bold", size = 12))
p <- p1 / p2
source("R/add_global_label.R")
p %>% add_global_label(Ylab ="% change in protected area", size=10, Ygap=0.04) 
ggsave("figures/Figure2_bins.png", dpi=1000, width=8, height=6)

sub_dat <- mar_dat %>% filter(path == "Temperature", marine == "coastal+marine", type=="Surface")
sub_dat_group <- mar_dat %>% filter(path == "Temperature", marine == "coastal+marine", type=="Surface") %>%
  group_by(var2, path, year_rcp) %>% summarise(perc_area=sum(perc_area)/n(), perc=sum(perc))
p1 <- ggplot() + geom_bar(data=sub_dat_group, aes(x=var2, y=perc_area*tot_area1e7), 
                    colour="black", fill="lightgrey", stat="identity") +
  geom_bar(data=sub_dat, aes(x = var2, y=perc*(perc_area/100)*tot_area1e7, fill=iucn_cat), 
           stat="identity", position = "stack") + 
  geom_bar(data=sub_dat_group, aes(x=var2, y=perc*(perc_area/100)*tot_area1e7), 
           col="black", fill="transparent", stat="identity") +
  geom_point(data=sub_dat_group, aes(x=var2, y=perc*tot_area1e7), col="red") +
  geom_line(data=sub_dat_group, aes(x=var2, y=perc*tot_area1e7), col="red", lty="dashed") + 
  facet_grid(year_rcp ~ path, scales="free_x", switch="x") + 
  labs(y="Area (1e07 km2)") + scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), 
                     sec.axis = sec_axis(~ ./tot_area1e7, breaks=c(0,2.5,5,7.5,10,12.5,15, 17.5))) + 
  scale_fill_manual(name="IUCN", values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"),
                    guide=guide_legend(reverse=T)) + 
  theme_bw() + theme(legend.position = "none", axis.title.x=element_blank(), 
                     strip.text = element_text(size=12, face="bold"), 
                     axis.title.y = element_text(size=12, face="bold"),
                     strip.text.y = element_blank(), strip.background = element_blank(), 
                     strip.placement="outside", panel.grid.minor = element_blank())

sub_dat <- mar_dat %>% filter(path == "Salinity", marine == "coastal+marine", type=="Surface")
sub_dat_group <- mar_dat %>% filter(path == "Salinity", marine == "coastal+marine", type=="Surface") %>%
  group_by(var2, path, year_rcp) %>% summarise(perc_area=sum(perc_area)/n(), perc=sum(perc))
p2 <- ggplot() + geom_bar(data=sub_dat_group, aes(x=var2, y=perc_area*tot_area1e7), 
                          colour="black", fill="lightgrey", stat="identity") +
  geom_bar(data=sub_dat, aes(x = var2, y=perc*(perc_area/100)*tot_area1e7, fill=iucn_cat), 
           stat="identity", position = "stack") + 
  geom_bar(data=sub_dat_group, aes(x=var2, y=perc*(perc_area/100)*tot_area1e7), 
           col="black", fill="transparent", stat="identity") +
  geom_point(data=sub_dat_group, aes(x=var2, y=perc*tot_area1e7), col="red") +
  geom_line(data=sub_dat_group, aes(x=var2, y=perc*tot_area1e7), col="red", lty="dashed") + 
  facet_grid(year_rcp ~ path, scales="free_x", switch="x") + 
  labs(y="Area (1e07 km2)") + scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), 
                     sec.axis = sec_axis(~ ./tot_area1e7, breaks=c(0,10,20,30,40))) + 
  scale_fill_manual(name="IUCN", values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"),
                    guide=guide_legend(reverse=T)) + 
  theme_bw() + theme(legend.position = "none", axis.title.x=element_blank(), 
                     strip.text = element_text(size=12, face="bold"), 
                     axis.title.y = element_blank(), strip.background = element_blank(), 
                     strip.placement="outside", panel.grid.minor = element_blank())

sub_dat <- mar_dat %>% filter(path == "Current.Velocity", marine == "coastal+marine", type=="Surface")
sub_dat_group <- mar_dat %>% filter(path == "Current.Velocity", marine == "coastal+marine", type=="Surface") %>% 
  group_by(var2, path, year_rcp) %>% summarise(perc_area=sum(perc_area)/n(), perc=sum(perc))
p3 <- ggplot() + geom_bar(data=sub_dat_group, aes(x=var2, y=perc_area*tot_area1e7), 
                          colour="black", fill="lightgrey", stat="identity") +
  geom_bar(data=sub_dat, aes(x = var2, y=perc*(perc_area/100)*tot_area1e7, fill=iucn_cat), 
           stat="identity", position = "stack") + 
  geom_bar(data=sub_dat_group, aes(x=var2, y=perc*(perc_area/100)*tot_area1e7), 
           col="black", fill="transparent", stat="identity") +
  geom_point(data=sub_dat_group, aes(x=var2, y=perc*tot_area1e7), col="red") +
  geom_line(data=sub_dat_group, aes(x=var2, y=perc*tot_area1e7), col="red", lty="dashed") + 
  facet_grid(year_rcp ~ path, scales="free_x", switch="x") + 
  labs(y="Area (1e07 km2)") + scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), 
                     sec.axis = sec_axis(~ ./tot_area1e7, breaks=c(0,5,10,15,20,25))) + 
  scale_fill_manual(name="IUCN", values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"),
                    guide=guide_legend(reverse=T)) + 
  theme_bw() + theme(legend.position = "none", axis.title.x=element_blank(), 
                     strip.text = element_text(size=12, face="bold"), 
                     axis.title.y = element_text(size=12, face="bold"),
                     strip.text.y = element_blank(), strip.background = element_blank(), 
                     strip.placement="outside", panel.grid.minor = element_blank())

sub_dat <- mar_dat %>% filter(path == "Ice.thickness", marine == "coastal+marine", type=="Surface")
sub_dat_group <- mar_dat %>% filter(path == "Ice.thickness", marine == "coastal+marine", type=="Surface") %>% 
  group_by(var2, path, year_rcp) %>% summarise(perc_area=sum(perc_area)/n(), perc=sum(perc))
p4 <- ggplot() + geom_bar(data=sub_dat_group, aes(x=var2, y=perc_area*tot_area1e7), 
                          colour="black", fill="lightgrey", stat="identity") +
  geom_bar(data=sub_dat, aes(x = var2, y=perc*(perc_area/100)*tot_area1e7, fill=iucn_cat), 
           stat="identity", position = "stack") + 
  geom_bar(data=sub_dat_group, aes(x=var2, y=perc*(perc_area/100)*tot_area1e7), 
           col="black", fill="transparent", stat="identity") +
  geom_point(data=sub_dat_group, aes(x=var2, y=perc*tot_area1e7), col="red") +
  geom_line(data=sub_dat_group, aes(x=var2, y=perc*tot_area1e7), col="red", lty="dashed") + 
  facet_grid(year_rcp ~ path, scales="free_x", switch="x") + 
  labs(y="") + scale_x_continuous(expand=c(0,0)) + 
  scale_y_continuous(expand=expansion(mult=c(0,.01)), 
                     sec.axis = sec_axis(~ ./tot_area1e7, breaks=c(0,20,40,60,80,100))) + 
  scale_fill_manual(name="IUCN", values=c("#D43F3AFF", "#EEA236FF", "#46B8DAFF", "#5CB85CFF"),
                    guide=guide_legend(reverse=T)) + 
  theme_bw() + theme(legend.position = "none", axis.title.x=element_blank(), 
                     strip.text = element_text(size=12, face="bold"), 
                     axis.title.y = element_text(size=12, face="bold"),
                     strip.placement="outside", strip.background = element_blank(), 
                     panel.grid.minor = element_blank())
(p1 | p2) / leg + plot_layout(heights=c(15,1))
ggsave("figures/pres_fut_hist_bins.png", dpi=1000, width=6, height=9)

(p3 | p4) / leg + plot_layout(heights=c(15,1)) #+ plot_annotation(tag_levels="a", tag_suffix=")")
ggsave("figures/pres_fut_hist_bins2.png", dpi=1000, width=6, height=9)
