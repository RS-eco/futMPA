#' ---
#' title: "Global Analysis of MPAs - Data summary"
#' author: "RS-eco"
#' ---

rm(list=ls()); gc()

# Set file directory
filedir <- "/home/matt/Documents/Bio-Oracle/allLayers"

#Automatically install required packages, which are not yet installed
packages <- c("raster", "SpaDES.tools", "tidyverse", "terra")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages); rm(new.packages)

# Load packages
l <- sapply(packages, require, character.only = TRUE, quietly=TRUE); rm(packages, l)

########################################

## Marine summary

# Stack mpa data to Bio-Oracle data
files <- list.files("extdata/", pattern="_bins.tif", full.names=T)
mpa_files <- list.files("extdata/", pattern="^mpa_cov.*\\.tif", full.names=T)

# stack() does not work as layers have different extents
# => Load raster layers as list
dat <- lapply(c(files, mpa_files), terra::rast)

# Find maximum extent of raster layers
ext <- lapply(dat, function(x) x %>% terra::ext() %>% as.vector() %>% as.data.frame)
ext <- bind_cols(ext)
ext <- t(ext) %>% as.data.frame()
colnames(ext) <- c("xmin", "xmax", "ymin", "ymax")
(ext_large <- ext %>% summarise(xmin=min(xmin), xmax=max(xmax), ymin=min(ymin), ymax=max(ymax)) %>% as.numeric() %>% terra::ext())

r <- terra::rast(resolution=c(8350, 10300), ext=terra::ext(ext_large), 
                    crs="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

# Extend and resample raster layers
dat <- lapply(dat, function(x){
  terra::expand(x, y=ext_large) %>% terra::resample(r) %>%
    terra::as.data.frame(xy=T)
  }); gc()
dat <- Reduce(function(...) dplyr::left_join(..., by=c("x","y"), all.x=TRUE), dat)
head(dat)
colnames(dat) <- sub(".Mean.BOv2_1_bins", "", sub(".Mean.Depth", "", sub(".Mean_bins", "", colnames(dat))))
saveRDS(dat, file="data/alldata_marine_df.rds", compress="xz")
rm(list=ls()); gc()
# Restart R session after this

dat <- readRDS("data/alldata_marine_df.rds")
colnames(dat)
head(dat)

# Select only variable colnames (without x,y and mpa values)
c_names <- colnames(dat)[3:57]

mar_dat_list <- lapply(c_names, function(x){
  # Summarise data and save to file
  sum_dat <- dat %>% dplyr::select(c(mpa_cov_I.II, mpa_cov_III.IV, mpa_cov_Notdesignated, mpa_cov_Total, mpa_cov_V.VI), matches(x)) %>%
    mutate_at(c("mpa_cov_I.II", "mpa_cov_III.IV", "mpa_cov_Notdesignated", "mpa_cov_Total", "mpa_cov_V.VI"), list(~ replace_na(., 0)))
  colnames(sum_dat)[ncol(sum_dat)] <- "var"
  sum_dat$var <- round(sum_dat$var, digits=0)
  sum_dat$path <- x
  head(sum_dat)
  mar_dat <- sum_dat %>% group_by(path, var) %>% summarise_all(list(~mean(., na.rm=T)))
  n_dat <- sum_dat %>% group_by(path, var) %>% summarise(n=n())
  mar_dat <- full_join(mar_dat, n_dat)
  return(mar_dat)
})
mar_dat <- bind_rows(mar_dat_list); rm(mar_dat_list); gc()
head(mar_dat)
tail(mar_dat)
mar_dat <- mar_dat %>% tidyr::separate(path, into=c("year", "rcp_type_path"), sep="[.]", extra="merge") 
mar_dat_fut2 <- mar_dat %>% filter(year == "X2100AOGCM") %>% tidyr::separate(rcp_type_path, into=c("rcp", "type", "path"), sep="[.]", extra="merge", fill="left")
mar_dat_fut <- mar_dat %>% filter(year == "X2050AOGCM") %>% tidyr::separate(rcp_type_path, into=c("rcp", "type", "path"), sep="[.]", extra="merge", fill="left")
mar_dat_pres <- mar_dat %>% filter(year == "Present") %>% tidyr::separate(rcp_type_path, into=c("type", "path"), sep="[.]", extra="merge", fill="left")
unique(mar_dat_pres$path)
mar_dat <- bind_rows(mar_dat_fut2, mar_dat_fut, mar_dat_pres)
head(mar_dat)
tail(mar_dat)
saveRDS(mar_dat, "data/summary_ind_marine_bins.rds", compress="xz")

dat <- dat %>% dplyr::select(c(mpa_cov_I.II, mpa_cov_III.IV, mpa_cov_Notdesignated, mpa_cov_Total, mpa_cov_V.VI), contains("Temperature"), contains("Salinity")); gc()
mar_dat2_list <- lapply(c("X2050AOGCM.RCP26.Benthic", "X2050AOGCM.RCP60.Benthic", "Present.Benthic", 
                          "X2050AOGCM.RCP26.Surface", "X2050AOGCM.RCP60.Surface", "Present.Surface",
                          "X2100AOGCM.RCP26.Benthic", "X2100AOGCM.RCP60.Benthic",
                          "X2100AOGCM.RCP26.Surface", "X2100AOGCM.RCP60.Surface"), function(x){
  sub_dat <- dat %>% mutate_at(c("mpa_cov_I.II", "mpa_cov_III.IV", "mpa_cov_Notdesignated", "mpa_cov_Total", "mpa_cov_V.VI"), list(~ replace_na(., 0))) %>%
    dplyr::select(c(mpa_cov_I.II, mpa_cov_III.IV, mpa_cov_Notdesignated, mpa_cov_Total, mpa_cov_V.VI), contains(x))
  mar_dat2 <- sub_dat %>% group_by_at(vars(-c(mpa_cov_I.II, mpa_cov_III.IV, mpa_cov_Notdesignated, mpa_cov_Total, mpa_cov_V.VI))) %>% 
    summarise_all(list(~mean(., na.rm=T)))
  n_dat2 <- sub_dat %>% group_by_at(vars(-c(mpa_cov_I.II, mpa_cov_III.IV, mpa_cov_Notdesignated, mpa_cov_Total, mpa_cov_V.VI))) %>% 
    summarise(n=n())
  mar_dat2 <- full_join(mar_dat2, n_dat2)
  colnames(mar_dat2)[1:2] <- c("temperature", "salinity")
  mar_dat2$temperature <- round(mar_dat2$temperature)
  mar_dat2$salinity <- round(mar_dat2$salinity)
  mar_dat2$type <- x
  return(mar_dat2)
})
mar_dat2 <- bind_rows(mar_dat2_list); rm(mar_dat2_list); gc()
head(mar_dat2)
mar_dat2 <- mar_dat2 %>% tidyr::separate(type, into=c("year", "rcp_type"), sep="[.]", extra="merge") %>%
  tidyr::separate(rcp_type, into=c("rcp", "type"), sep="[.]", extra="merge", fill="left")
head(mar_dat2)
tail(mar_dat2)
saveRDS(mar_dat2, "data/summary_all_marine_bins.rds", compress="xz")

########################################

## Coastal summary

# Stack cpa data to Bio-Oracle data
files <- list.files("extdata/", pattern="_bins.tif", full.names=T)
cpa_files <- list.files("extdata/", pattern="^cpa_cov.*\\.tif", full.names=T)

# Does not work as layers have different extents

# Load raster layers as list
dat <- lapply(c(files, cpa_files), raster::raster)

# Find maximum extent of raster layers
ext <- lapply(dat, function(x) x %>% extent %>% as.vector() %>% as.data.frame)
ext <- bind_cols(ext)
ext <- t(ext) %>% as.data.frame()
colnames(ext) <- c("xmin", "xmax", "ymin", "ymax")
(ext_large <- ext %>% summarise(xmin=min(xmin), xmax=max(xmax), ymin=min(ymin), ymax=max(ymax)) %>% as.numeric())

r <- raster::raster(resolution=c(8350, 10300), ext=extent(ext_large), 
                    crs="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

# Extend and resample raster layers
dat <- lapply(dat, function(x) extend(x, y=ext_large))
dat <- lapply(dat, function(x){raster::resample(x, r)})
dat <- raster::stack(dat); gc()
dat <- as.data.frame(raster::rasterToPoints(dat)); gc()
colnames(dat) <- sub(".Mean.BOv2_1_bins", "", sub(".Mean.Depth", "", sub(".Mean_bins", "", colnames(dat))))
saveRDS(dat, file="data/alldata_coastal_df.rds", compress="xz")
rm(list=ls()); gc()
# Restart R session after this

dat <- readRDS("data/alldata_coastal_df.rds")
colnames(dat)

# Select only variable colnames (without x,y and mpa values)
(c_names <- colnames(dat)[3:57])

mar_dat_list <- lapply(c_names, function(x){
  
  # Summarise data and save to file
  sum_dat <- dat %>% mutate_at(c("cpa_cov_I.II", "cpa_cov_III.IV", "cpa_cov_Notdesignated", "cpa_cov_Total", "cpa_cov_V.VI"), list(~ replace_na(., 0))) %>%
    dplyr::select(c(cpa_cov_I.II, cpa_cov_III.IV, cpa_cov_Notdesignated, cpa_cov_Total, cpa_cov_V.VI), matches(x))
  colnames(sum_dat)[ncol(sum_dat)] <- "var"
  sum_dat$var <- round(sum_dat$var, digits=0)
  sum_dat$path <- x
  head(sum_dat)
  mar_dat <- sum_dat %>% group_by(path, var) %>% summarise_all(list(~mean(., na.rm=T)))
  n_dat <- sum_dat %>% group_by(path, var) %>% summarise(n=n())
  mar_dat <- full_join(mar_dat, n_dat)
  return(mar_dat)
})
mar_dat <- bind_rows(mar_dat_list); rm(mar_dat_list); gc()
mar_dat <- mar_dat %>% tidyr::separate(path, into=c("year", "rcp_type_path"), sep="[.]", extra="merge") 
mar_dat_fut2 <- mar_dat %>% filter(year == "X2100AOGCM") %>% tidyr::separate(rcp_type_path, into=c("rcp", "type", "path"), sep="[.]", extra="merge", fill="left")
mar_dat_fut <- mar_dat %>% filter(year == "X2050AOGCM") %>% tidyr::separate(rcp_type_path, into=c("rcp", "type", "path"), sep="[.]", extra="merge", fill="left")
mar_dat_pres <- mar_dat %>% filter(year == "Present") %>% tidyr::separate(rcp_type_path, into=c("type", "path"), sep="[.]", extra="merge", fill="left")
unique(mar_dat_pres$path)
mar_dat <- bind_rows(mar_dat_fut2, mar_dat_fut, mar_dat_pres)
head(mar_dat)
tail(mar_dat)
saveRDS(mar_dat, "data/summary_ind_coast_bins.rds", compress="xz")

dat <- dat %>% dplyr::select(c(cpa_cov_I.II, cpa_cov_III.IV, cpa_cov_Notdesignated, cpa_cov_Total, cpa_cov_V.VI), contains("Temperature"), contains("Salinity")); gc()
mar_dat2_list <- lapply(c("X2050AOGCM.RCP26.Benthic", "X2050AOGCM.RCP60.Benthic", "Present.Benthic", 
                          "X2050AOGCM.RCP26.Surface", "X2050AOGCM.RCP60.Surface", "Present.Surface",
                          "X2100AOGCM.RCP26.Benthic", "X2100AOGCM.RCP60.Benthic",
                          "X2100AOGCM.RCP26.Surface", "X2100AOGCM.RCP60.Surface"), function(x){
  sub_dat <- dat %>% mutate_at(c("cpa_cov_I.II", "cpa_cov_III.IV", "cpa_cov_Notdesignated", "cpa_cov_Total", "cpa_cov_V.VI"), list(~ replace_na(., 0))) %>%
    dplyr::select(c(cpa_cov_I.II, cpa_cov_III.IV, cpa_cov_Notdesignated, cpa_cov_Total, cpa_cov_V.VI), contains(x))
  mar_dat2 <- sub_dat %>% group_by_at(vars(-c(cpa_cov_I.II, cpa_cov_III.IV, cpa_cov_Notdesignated, cpa_cov_Total, cpa_cov_V.VI))) %>% 
    summarise_all(list(~mean(., na.rm=T)))
  n_dat2 <- sub_dat %>% group_by_at(vars(-c(cpa_cov_I.II, cpa_cov_III.IV, cpa_cov_Notdesignated, cpa_cov_Total, cpa_cov_V.VI))) %>% 
    summarise(n=n())
  mar_dat2 <- full_join(mar_dat2, n_dat2)
  colnames(mar_dat2)[1:2] <- c("temperature", "salinity")
  mar_dat2$temperature <- round(mar_dat2$temperature)
  mar_dat2$salinity <- round(mar_dat2$salinity)
  mar_dat2$type <- x
  return(mar_dat2)
})
mar_dat2 <- bind_rows(mar_dat2_list); rm(mar_dat2_list); gc()
head(mar_dat2)
mar_dat2 <- mar_dat2 %>% tidyr::separate(type, into=c("year", "rcp_type"), sep="[.]", extra="merge") %>%
  tidyr::separate(rcp_type, into=c("rcp", "type"), sep="[.]", extra="merge", fill="left")
head(mar_dat2)
tail(mar_dat2)
saveRDS(mar_dat2, "data/summary_all_coast_bins.rds", compress="xz")

########################################

## Coastal & marine summary

# Stack cpa data to Bio-Oracle data
files <- list.files("extdata/", pattern="_bins.tif", full.names=T)
pa_files <- list.files("extdata/", pattern="^pa_cov.*\\.tif", full.names=T)

# Does not work as layers have different extents

# Load raster layers as list
dat <- lapply(c(files, pa_files), raster::raster)

# Find maximum extent of raster layers
ext <- lapply(dat, function(x) x %>% extent %>% as.vector() %>% as.data.frame)
ext <- bind_cols(ext)
ext <- t(ext) %>% as.data.frame()
colnames(ext) <- c("xmin", "xmax", "ymin", "ymax")
(ext_large <- ext %>% summarise(xmin=min(xmin), xmax=max(xmax), ymin=min(ymin), ymax=max(ymax)) %>% as.numeric())

r <- raster::raster(resolution=c(8350, 10300), ext=extent(ext_large), 
                    crs="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

# Extend and resample raster layers
dat <- lapply(dat, function(x) extend(x, y=ext_large)); gc()
dat <- lapply(dat, function(x){raster::resample(x, r)}); gc()
dat <- raster::stack(dat); gc()
dat <- as.data.frame(raster::rasterToPoints(dat)); gc()
colnames(dat) <- sub(".Mean.BOv2_1_bins", "", sub(".Mean.Depth", "", sub(".Mean_bins", "", colnames(dat))))
saveRDS(dat, file="data/alldata_coastal+marine_df.rds", compress="xz")
rm(list=ls()); gc()
# Restart R session after this

dat <- readRDS("data/alldata_coastal+marine_df.rds")
colnames(dat)

# Select only variable colnames (without x,y and mpa values)
(c_names <- colnames(dat)[3:57])

mar_dat_list <- lapply(c_names, function(x){
  # Summarise data and save to file
  sum_dat <- dat %>% mutate_at(c("pa_cov_I.II", "pa_cov_III.IV", "pa_cov_Notdesignated", "pa_cov_Total", "pa_cov_V.VI"), list(~ replace_na(., 0))) %>%
    dplyr::select(c(pa_cov_I.II, pa_cov_III.IV, pa_cov_Notdesignated, pa_cov_Total, pa_cov_V.VI), matches(x))
  colnames(sum_dat)[ncol(sum_dat)] <- "var"
  sum_dat$var <- round(sum_dat$var, digits=0)
  sum_dat$path <- x
  head(sum_dat)
  mar_dat <- sum_dat %>% group_by(path, var) %>% summarise_all(list(~mean(., na.rm=T)))
  n_dat <- sum_dat %>% group_by(path, var) %>% summarise(n=n())
  mar_dat <- full_join(mar_dat, n_dat)
  return(mar_dat)
})
mar_dat <- bind_rows(mar_dat_list); rm(mar_dat_list); gc()
mar_dat <- mar_dat %>% tidyr::separate(path, into=c("year", "rcp_type_path"), sep="[.]", extra="merge") 
mar_dat_fut2 <- mar_dat %>% filter(year == "X2100AOGCM") %>% tidyr::separate(rcp_type_path, into=c("rcp", "type", "path"), sep="[.]", extra="merge", fill="left")
mar_dat_fut <- mar_dat %>% filter(year == "X2050AOGCM") %>% tidyr::separate(rcp_type_path, into=c("rcp", "type", "path"), sep="[.]", extra="merge", fill="left")
mar_dat_pres <- mar_dat %>% filter(year == "Present") %>% tidyr::separate(rcp_type_path, into=c("type", "path"), sep="[.]", extra="merge", fill="left")
unique(mar_dat_pres$path)
mar_dat <- bind_rows(mar_dat_fut2, mar_dat_fut, mar_dat_pres)
head(mar_dat)
tail(mar_dat)
saveRDS(mar_dat, "data/summary_ind_coastal+marine_bins.rds", compress="xz")

dat <- dat %>% dplyr::select(c(pa_cov_I.II, pa_cov_III.IV, pa_cov_Notdesignated, pa_cov_Total, pa_cov_V.VI), contains("Temperature"), contains("Salinity")); gc()
mar_dat2_list <- lapply(c("X2050AOGCM.RCP26.Benthic", "X2050AOGCM.RCP60.Benthic", "Present.Benthic", 
                          "X2050AOGCM.RCP26.Surface", "X2050AOGCM.RCP60.Surface", "Present.Surface",
                          "X2100AOGCM.RCP26.Benthic", "X2100AOGCM.RCP60.Benthic",
                          "X2100AOGCM.RCP26.Surface", "X2100AOGCM.RCP60.Surface"), function(x){
                            sub_dat <- dat %>% mutate_at(c("pa_cov_I.II", "pa_cov_III.IV", "pa_cov_Notdesignated", "pa_cov_Total", "pa_cov_V.VI"), list(~ replace_na(., 0))) %>%
                              dplyr::select(c(pa_cov_I.II, pa_cov_III.IV, pa_cov_Notdesignated, pa_cov_Total, pa_cov_V.VI), contains(x))
                            mar_dat2 <- sub_dat %>% group_by_at(vars(-c(pa_cov_I.II, pa_cov_III.IV, pa_cov_Notdesignated, pa_cov_Total, pa_cov_V.VI))) %>% 
                              summarise_all(list(~mean(., na.rm=T)))
                            n_dat2 <- sub_dat %>% group_by_at(vars(-c(pa_cov_I.II, pa_cov_III.IV, pa_cov_Notdesignated, pa_cov_Total, pa_cov_V.VI))) %>% 
                              summarise(n=n())
                            mar_dat2 <- full_join(mar_dat2, n_dat2)
                            colnames(mar_dat2)[1:2] <- c("temperature", "salinity")
                            mar_dat2$temperature <- round(mar_dat2$temperature)
                            mar_dat2$salinity <- round(mar_dat2$salinity)
                            mar_dat2$type <- x
                            return(mar_dat2)
                          })
mar_dat2 <- bind_rows(mar_dat2_list); rm(mar_dat2_list); gc()
head(mar_dat2)
mar_dat2 <- mar_dat2 %>% tidyr::separate(type, into=c("year", "rcp_type"), sep="[.]", extra="merge") %>%
  tidyr::separate(rcp_type, into=c("rcp", "type"), sep="[.]", extra="merge", fill="left")
head(mar_dat2)
tail(mar_dat2)
saveRDS(mar_dat2, "data/summary_all_coastal+marine_bins.rds", compress="xz")

########################################