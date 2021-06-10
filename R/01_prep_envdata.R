
rm(list=ls()); gc()

# Set file directory
filedir <- "/home/matt/Documents/Bio-Oracle/allLayers"

# Load packages
library(tidyverse); library(terra)

# List available files
files <- list.files(filedir, pattern="Mean.tif", full.names=T)

####################

# One file (calcite) is in lower letters
lower_files <- list.files(filedir, pattern="mean.tif", full.names=T)
files <- c(files, lower_files)

####################

# Define combinations of layers and time_periods
layers <- c("Surface", "Benthic.Mean.Depth")
time_periods <- c("Present", "2050AOGCM.RCP26", "2050AOGCM.RCP60", "2100AOGCM.RCP26", "2100AOGCM.RCP60")

combinations <- expand.grid(time_periods, layers) %>% unite("Var", Var1:Var2, sep=".") %>% unlist()

# Identify variables per time_period and layer
vars <- lapply(combinations, function(x){
  sapply(sub(paste0(x, "."), "", basename(files[grepl(files, pattern=x)])), function(y){sub(".Mean", "", paste0(strsplit(y, split="[.]")[[1]][1:2], collapse="."))})
})
vars

####################

# Transform data into CEA projection

lapply(1:length(files), function(z){
  if(!file.exists(paste0("extdata/", gsub(".tif", "", basename(files[z])), "_wm.tif"))){
    dat <- raster::raster(files[z])
    dat
    if(is.na(raster::projection(dat))){
      raster::projection(dat) <- "+proj=longlat +datum=WGS84 +no_defs"
      print(paste0(z, "is missing projection"))
    }
    if(raster::projection(dat) != "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"){
      raster::projectRaster(dat, crs="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs", 
                            filename=paste0("extdata/", gsub(".tif", "", basename(files[z])), "_wm.tif"), 
                            format="GTiff", overwrite=T)
      return(NULL)
    } else{
      #writeRaster(dat, filename=paste0("extdata/", gsub(".tif", "", basename(files[z])), "_cea.tif"), 
      #            format="GTiff", overwrite=T)
      print(paste0(z, "has already CEA projection"))
    }
  }
})

#' **Note:** Some files are in cea rather than wgs84 projection and do not fit with the data extent of the other variables, 
#' hence we leave them out here.

####################

# Identify range of individual data files

if(!file.exists("data/min_max_files.rds")){
  files <- list.files("extdata/", pattern="_wm.tif", full.names=T)
  sum_dat <- lapply(1:length(files), function(z){
    dat <- terra::rast(files[z])
    df <- data.frame(filename=files[z], min=terra::global(dat, "min", na.rm=T), 
                     max=terra::global(dat, "max", na.rm=T))
    return(df)
  })
  sum_dat <- bind_rows(sum_dat)
  rownames(sum_dat) <- NULL
  saveRDS(sum_dat, file="data/min_max_files.rds", compress="xz")
} else{
  sum_dat <- readRDS("data/min_max_files.rds")
}

# Load data and define variables
sum_dat <- sum_dat %>% 
  rowwise() %>% 
  mutate(var = sub(".tif", "", 
                   sub("_wm", "", 
                       sub(".Mean", "", 
                           sub("Benthic.", "", 
                               sub("Surface.", "", 
                                   paste0(strsplit(sub(".Mean.Depth", "", 
                                                       basename(filename)), 
                                                   split="[.]")[[1]][3:5],
                                          collapse=".")))))))
head(sum_dat, 7)

# Identify unique variables
unique(sum_dat$var)

# Create overall summary per variable
sum_var <- sum_dat %>% group_by(var) %>%
  summarise(min=floor(min(min)), max=ceiling(max(max)))

####################

# Split data into equally-sized bins

sum_var

# Create matrix for relassfication
m <- c(list(as.matrix(data.frame(x=seq(0, 6.75, by=0.25), y=seq(0.25, 7, by=0.25), z=1:28)),
            as.matrix(data.frame(x=seq(-5, 1.75, by=0.25), y=seq(-4.75, 2, by=0.25), z=1:28)),
            as.matrix(data.frame(x=seq(0, 1.9, by=0.1), y=seq(0.1, 2, by=0.1), z=1:20)),
            as.matrix(data.frame(x=seq(-1, 17.5, by=0.5), y=seq(-0.5, 18, by=0.5), z=1:38)),
            as.matrix(data.frame(x=seq(0, 413, by=7), y=seq(7, 420, by=7), z=1:60)),
            as.matrix(data.frame(x=seq(-1, 0.9, by=0.1), y=seq(-0.9, 1, by=0.1), z=1:20)),
            as.matrix(data.frame(x=seq(-1, 7.5, by=0.5), y=seq(-0.5, 8, by=0.5), z=1:18)),
            as.matrix(data.frame(x=seq(0, 0.95, by=0.05), y=seq(0.05, 1, by=0.05), z=1:20)),
            as.matrix(data.frame(x=seq(-1, 53, by=1), y=seq(0, 54, by=1), z=1:55)),
            as.matrix(data.frame(x=seq(-1, 107, by=2), y=seq(1, 109, by=2), z=1:55)),
            as.matrix(data.frame(x=seq(0, 3.8, by=0.2), y=seq(0.2, 4, by=0.2), z=1:20)),
            as.matrix(data.frame(x=seq(0, 19, by=1), y=seq(1, 20, by=1), z=1:20)),
            as.matrix(data.frame(x=seq(-1, 0.9, by=0.1), y=seq(-0.9, 1, by=0.1), z=1:20)),
            as.matrix(data.frame(x=seq(0, 41, by=1), y=seq(1, 42, by=1), z=1:42)),
            as.matrix(data.frame(x=seq(0, 310, by=10), y=seq(10, 320, by=10), z=1:32)),
            as.matrix(data.frame(x=seq(-2, 33.5, by=0.5), y=seq(-1.5, 34, by=0.5), z=1:72))))
names(m) <- unique(sum_var$var)

# Reclassify data layers

lapply(1:nrow(sum_dat), function(x){
  dat <- terra::rast(sum_dat$filename[x])
  rc2 <- terra::classify(dat, rcl=m[[sum_dat$var[x]]],
                         filename=paste0("extdata/", sub("_wm.tif", "_bins.tif", basename(sum_dat$filename[x]))), 
                         format="GTiff", overwrite=T)
})

dat <- terra::rast(paste0("extdata/", sub("_wm.tif", "_bins.tif", basename(sum_dat$filename[1]))))
dat

####################

# Split data into equal-frequency bins

### Summarise raster stack
if(!file.exists("data/summary_mar_perc.rds")){
  files <- list.files("extdata/", pattern="_wm.tif", full.names=T)
  summary_mar <- lapply(files, function(k){
    dat <- terra::rast(k)
    x <- terra::values(dat)
    x <- na.omit(x); gc()
    bins_perc <- binr::bins(x, target.bins=25, max.breaks=25)
    bins_perc_val <- binr::bins.getvals(bins_perc)
    rclmat <- data.frame(var=sub("_wm.tif", "", basename(k)),
                         x=as.numeric(c(attr(bins_perc_val, "binlo")[1:(length(bins_perc$binct))])), 
                         y=as.numeric(c(attr(bins_perc_val, "binlo")[2:(length(bins_perc$binct))],attr(bins_perc_val, "binhi")[length(bins_perc$binct)])), 
                         z=1:length(bins_perc$binct))
    rclmat[nrow(rclmat),3] <- rclmat[nrow(rclmat),3]+1
    return(rclmat)
  })
  summary_mar <- dplyr::bind_rows(summary_mar)
  saveRDS(summary_mar, file="data/summary_mar_perc.rds", compress="xz")
  gc() 
} else{
  summary_mar <- readRDS("data/summary_mar_perc.rds")
}


# Reclassify

#files <- c("extdata/bathy_30s_wm.tif", list.files("extdata", pattern="biogeo.*\\_wm.tif", full.names=T))
#m_ee <- readRDS("data/summary_marspec_perc_optim.rds")
#m_ee[1,2] <- m_ee[1,2]-1
#m_ee[185,2] <- m_ee[185,2]-1
#m_ee[335,2] <- m_ee[335,2]-1
#m_ee[431,2] <- floor(m_ee[431,2])
#m_ee[528,2] <- floor(m_ee[528,2])
#m <- m_ee %>% group_split(var)

# Reclassify data
#lapply(1:length(files), function(x){
#  dat <- terra::rast(files[x])
#  rc2 <- terra::classify(dat, rcl=as.matrix(m[[x]][,2:4]), right=T, include.lowest=T, 
#                         filename=paste0("extdata/", sub("_wm.tif", "_perc.tif", basename(files[x]))), 
#                         format="GTiff", overwrite=T); rm(dat); gc()
#})

#dat <- terra::rast(paste0("extdata/", sub("_wm.tif", "_perc.tif", basename(files[7]))))
#dat

####################