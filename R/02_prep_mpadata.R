#' ---
#' title: "Global Analysis of Protected Areas - Data preparation"
#' author: "RS-eco"
#' ---

# Set working directory
workdir <- "/home/matt/Documents/futMPA"
setwd(workdir)

# Set file directory
filedir <- "/home/matt/Documents/"

#Automatically install required packages, which are not yet installed
packages <- c("raster", "fasterize", "sf", "tidyverse", "SpaDES.tools")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="https://cloud.r-project.org"); rm(new.packages)

# Load packages
l <- sapply(packages, require, character.only = TRUE, quietly=TRUE); rm(packages, l)

# Define IUCN category split
iucn_cat2 <- c("I-II", "III-IV", "V-VI", "Not designated", "Total")

########################################

#' ## Get WDPA Data and subset by IUCN category for marine and terrestrial separately

#' Download Geodatabase of all PAs
#download.file("http://d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_wdpa_shp.zip",
#              destfile=paste0(filedir, "WDPA/WDPA_WDOECM_wdpa_shp.zip"))

# Downloaded 9th of December 2020

#' Unzip shapefile
#unzip(paste0(filedir, "WDPA/WDPA_WDOECM_wdpa_shp.zip"), exdir=paste0(filedir, "WDPA"))

#' Divide protected areas in marine and terrestrial protected areas, always keep Coastal
#' and save new sf object as RDS

wdpa <- sf::st_read(dsn=paste0(filedir, "WDPA/WDPA_WDOECM_wdpa_shp0/WDPA_WDOECM_wdpa_shp-polygons.shp"))
wdpa1 <- sf::st_read(dsn=paste0(filedir, "WDPA/WDPA_WDOECM_wdpa_shp1/WDPA_WDOECM_wdpa_shp-polygons.shp"))
wdpa2 <- sf::st_read(dsn=paste0(filedir, "WDPA/WDPA_WDOECM_wdpa_shp2/WDPA_WDOECM_wdpa_shp-polygons.shp"))
wdpa <- bind_rows(wdpa, wdpa1, wdpa2); rm(wdpa1, wdpa2); gc()

wdpa_point <- sf::st_read(dsn=paste0(filedir, "WDPA/WDPA_WDOECM_wdpa_shp0/WDPA_WDOECM_wdpa_shp-points.shp"))
wdpa_point1 <- sf::st_read(dsn=paste0(filedir, "WDPA/WDPA_WDOECM_wdpa_shp1/WDPA_WDOECM_wdpa_shp-points.shp"))
wdpa_point2 <- sf::st_read(dsn=paste0(filedir, "WDPA/WDPA_WDOECM_wdpa_shp2/WDPA_WDOECM_wdpa_shp-points.shp"))
wdpa_point <- bind_rows(wdpa_point, wdpa_point1, wdpa_point2); rm(wdpa_point1, wdpa_point2); gc()

# Marine
mpa <- wdpa[wdpa$MARINE == 2,]
saveRDS(mpa, "extdata/WDPA_Dec2020_Marine_Total.rds", compress="xz")

mpa_point <- wdpa_point[wdpa_point$MARINE == 2,]
saveRDS(mpa_point, "extdata/WDPA_Dec2020_Marine_Point.rds", compress="xz")

# Coastal
cpa <- wdpa[wdpa$MARINE == 1,]
saveRDS(cpa, "extdata/WDPA_Dec2020_Coastal_Total.rds", compress="xz")

cpa_point <- wdpa_point[wdpa_point$MARINE == 1,]
saveRDS(cpa_point, "extdata/WDPA_Dec2020_Coastal_Point.rds", compress="xz")

# Coastal & Marine
pa <- wdpa[wdpa$MARINE %in% c(1,2),]
saveRDS(pa, "extdata/WDPA_Dec2020_All_Total.rds", compress="xz")

pa_point <- wdpa_point[wdpa_point$MARINE %in% c(1,2),]
saveRDS(pa_point, "extdata/WDPA_Dec2020_All_Point.rds", compress="xz")

rm(list=ls()); gc()
#' Restart R after this

# Process marine data
mpa <- readRDS("extdata/WDPA_Dec2020_Marine_Total.rds")
mpa_point <- readRDS("extdata/WDPA_Dec2020_Marine_Point.rds")
lapply(iucn_cat2[1:4], function(k){
  if(!file.exists(paste0("extdata/WDPA_Dec2020_Marine_", sub(" ", "", k), ".rds"))){
    if(k == "I-II"){
      mpa_sub <- mpa[mpa$IUCN_CAT %in% c("Ia", "Ib", "II"),]; gc()
      mpa_point_sub <- mpa_point[mpa_point$IUCN_CAT %in% c("Ia", "Ib", "II"),]
    } else if(k == "III-IV"){
      mpa_sub <- mpa[mpa$IUCN_CAT %in% c("III", "IV"),]; gc()
      mpa_point_sub <- mpa_point[mpa_point$IUCN_CAT %in% c("III", "IV"),]
    } else if(k == "V-VI"){
      mpa_sub <- mpa[mpa$IUCN_CAT %in% c("V","VI"),]; gc()
      mpa_point_sub <- mpa_point[mpa_point$IUCN_CAT %in% c("V", "VI"),]
    } else if(k == "Not designated"){
      mpa_sub <- mpa[mpa$IUCN_CAT %in% c("Not Reported", "Not Applicable", "Not Assigned"),]; gc()
      mpa_point_sub <- mpa_point[mpa_point$IUCN_CAT %in% c("Not Reported", "Not Applicable", "Not Assigned"),]
    }
    saveRDS(mpa_sub, paste0("extdata/WDPA_Dec2020_Marine_", sub(" ", "", k), ".rds"), compress="xz"); rm(mpa_sub)
    saveRDS(mpa_point_sub, paste0("extdata/WDPA_Dec2020_Marine_point_", sub(" ", "", k), ".rds"), compress="xz"); rm(mpa_point_sub)
    invisible(gc())
  }
}); rm(mpa, mpa_point); gc()

# Process coastal data
cpa <- readRDS("extdata/WDPA_Dec2020_Coastal_Total.rds")
cpa_point <- readRDS("extdata/WDPA_Dec2020_Coastal_Point.rds")
lapply(iucn_cat2[1:4], function(k){
  if(!file.exists(paste0("extdata/WDPA_Dec2020_Coastal_", sub(" ", "", k), ".rds"))){
    if(k == "I-II"){
      cpa_sub <- cpa[cpa$IUCN_CAT %in% c("Ia", "Ib", "II"),]; gc()
      cpa_point_sub <- cpa_point[cpa_point$IUCN_CAT %in% c("Ia", "Ib", "II"),]
    } else if(k == "III-IV"){
      cpa_sub <- cpa[cpa$IUCN_CAT %in% c("III", "IV"),]; gc()
      cpa_point_sub <- cpa_point[cpa_point$IUCN_CAT %in% c("III", "IV"),]
    } else if(k == "V-VI"){
      cpa_sub <- cpa[cpa$IUCN_CAT %in% c("V","VI"),]; gc()
      cpa_point_sub <- cpa_point[cpa_point$IUCN_CAT %in% c("V", "VI"),]
    } else if(k == "Not designated"){
      cpa_sub <- cpa[cpa$IUCN_CAT %in% c("Not Reported", "Not Applicable", "Not Assigned"),]; gc()
      cpa_point_sub <- cpa_point[cpa_point$IUCN_CAT %in% c("Not Reported", "Not Applicable", "Not Assigned"),]
    }
    saveRDS(cpa_sub, paste0("extdata/WDPA_Dec2020_Coastal_", sub(" ", "", k), ".rds"), compress="xz"); rm(cpa_sub)
    saveRDS(cpa_point_sub, paste0("extdata/WDPA_Dec2020_Coastal_point_", sub(" ", "", k), ".rds"), compress="xz"); rm(cpa_point_sub)
    invisible(gc())
  }
}); rm(cpa, cpa_point); gc()

# Process coastal & marine data
pa <- readRDS("extdata/WDPA_Dec2020_All_Total.rds")
pa_point <- readRDS("extdata/WDPA_Dec2020_All_Point.rds")
lapply(iucn_cat2[1:4], function(k){
  if(!file.exists(paste0("extdata/WDPA_Dec2020_All_", sub(" ", "", k), ".rds"))){
    if(k == "I-II"){
      pa_sub <- pa[pa$IUCN_CAT %in% c("Ia", "Ib", "II"),]; gc()
      pa_point_sub <- pa_point[pa_point$IUCN_CAT %in% c("Ia", "Ib", "II"),]
    } else if(k == "III-IV"){
      pa_sub <- pa[pa$IUCN_CAT %in% c("III", "IV"),]; gc()
      pa_point_sub <- pa_point[pa_point$IUCN_CAT %in% c("III", "IV"),]
    } else if(k == "V-VI"){
      pa_sub <- pa[pa$IUCN_CAT %in% c("V","VI"),]; gc()
      pa_point_sub <- pa_point[pa_point$IUCN_CAT %in% c("V", "VI"),]
    } else if(k == "Not designated"){
      pa_sub <- pa[pa$IUCN_CAT %in% c("Not Reported", "Not Applicable", "Not Assigned"),]; gc()
      pa_point_sub <- pa_point[pa_point$IUCN_CAT %in% c("Not Reported", "Not Applicable", "Not Assigned"),]
    }
    saveRDS(pa_sub, paste0("extdata/WDPA_Dec2020_All_", sub(" ", "", k), ".rds"), compress="xz"); rm(pa_sub)
    saveRDS(pa_point_sub, paste0("extdata/WDPA_Dec2020_All_point_", sub(" ", "", k), ".rds"), compress="xz"); rm(pa_point_sub)
    invisible(gc())
  }
}); rm(pa, pa_point); gc()

########################################

#' ## Rasterize data according to area covered by each IUCN category

# Get coverage of mpa data
lapply(iucn_cat2, function(k){
  if(!file.exists(paste0("extdata/mpa_cov_", sub(" ", "", k), ".tif"))){
    if(file.exists(paste0("extdata/WDPA_Dec2020_Marine_", sub(" ", "", k), ".rds"))){
      mpa <- readRDS(paste0("extdata/WDPA_Dec2020_Marine_", sub(" ", "", k), ".rds"))
      
      # Specify raster resolution (1km)
      r <- raster::raster(nrow=21600, ncol=43200, xmn=-180, xmx=180, ymn=-90, ymx=90,
                          crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      r_sub <- raster::crop(r, as(mpa, "Spatial"))
      tmpdir <- file.path(tempdir(), "splitRaster")
      dir.create(tmpdir)
      r_sub <- SpaDES.tools::splitRaster(r_sub, nx=5, ny=5, path=file.path(tmpdir, "mpa_sub"))
      lapply(1:length(r_sub), function(l){
        if(!file.exists(paste0("data/mpa_cov_", sub(" ", "", k), "_", l, ".tif"))){
          mpa_cov <- fasterize::fasterize(mpa, r_sub[[l]])
          mpa_cov <- raster::aggregate(mpa_cov, fact=10, fun=sum, 
                                       filename=paste0("extdata/mpa_cov_", sub(" ", "", k), "_", l, ".tif"), 
                                       format="GTiff", options="COMPRESS=LZW", overwrite=T) 
        }
      }); rm(mpa, r_sub); gc()
    }
  }
})

lapply(iucn_cat2, function(k){
  if(!file.exists(paste0("extdata/mpa_cov_", sub(" ", "", k), ".tif"))){
    mpa_files <- list.files(path="extdata", pattern=paste0("mpa_cov_", sub(" ", "", k), "_"), full.names=T)
    if(length(mpa_files) == 25){
      mpa_cov <- lapply(mpa_files, raster)
      mpa_cov <- do.call(raster::merge, c(mpa_cov, tolerance=0.5))
      r <- raster::raster(nrow=2160, ncol=4320, xmn=-180, xmx=180, ymn=-90, ymx=90,
                          crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      mpa_cov <- raster::extend(mpa_cov, r)
      mpa_cov <- raster::projectRaster(mpa_cov, crs="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs", 
                                       filename=paste0("extdata/mpa_cov_", sub(" ", "", k), ".tif"), 
                                       format="GTiff", options="COMPRESS=LZW", overwrite=T); file.remove(mpa_files); rm(mpa_cov); gc()
    }
  }
})

dat <- raster::raster("extdata/mpa_cov_I-II.tif")
dat
rm(list=ls());gc()

# Get coverage of cpa data
lapply(iucn_cat2, function(k){
  if(!file.exists(paste0("extdata/cpa_cov_", sub(" ", "", k), ".tif"))){
    if(file.exists(paste0("extdata/WDPA_Dec2020_Coastal_", sub(" ", "", k), ".rds"))){
      cpa <- readRDS(paste0("extdata/WDPA_Dec2020_Coastal_", sub(" ", "", k), ".rds"))
      
      # Specify raster resolution (1km)
      r <- raster::raster(nrow=21600, ncol=43200, xmn=-180, xmx=180, ymn=-90, ymx=90,
                          crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      r_sub <- raster::crop(r, as(cpa, "Spatial"))
      tmpdir <- file.path(tempdir(), "splitRaster")
      dir.create(tmpdir)
      r_sub <- SpaDES.tools::splitRaster(r_sub, nx=5, ny=5, path=file.path(tmpdir, "mpa_sub"))
      lapply(1:length(r_sub), function(l){
        if(!file.exists(paste0("data/cpa_cov_", sub(" ", "", k), "_", l, ".tif"))){
          cpa_cov <- fasterize::fasterize(cpa, r_sub[[l]])
          cpa_cov <- raster::aggregate(cpa_cov, fact=10, fun=sum, 
                                       filename=paste0("extdata/cpa_cov_", sub(" ", "", k), "_", l, ".tif"), 
                                       format="GTiff", options="COMPRESS=LZW", overwrite=T) 
        }
      }); rm(cpa, r_sub); gc()
    }
  }
})

lapply(iucn_cat2, function(k){
  if(!file.exists(paste0("extdata/cpa_cov_", sub(" ", "", k), ".tif"))){
    cpa_files <- list.files(path="extdata", pattern=paste0("cpa_cov_", sub(" ", "", k), "_"), full.names=T)
    if(length(cpa_files) == 25){
      cpa_cov <- lapply(cpa_files, raster)
      cpa_cov <- do.call(raster::merge, c(cpa_cov, tolerance=0.5))
      r <- raster::raster(nrow=2160, ncol=4320, xmn=-180, xmx=180, ymn=-90, ymx=90,
                          crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      cpa_cov <- raster::extend(cpa_cov, r)
      cpa_cov <- raster::projectRaster(cpa_cov, crs="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs", 
                                       filename=paste0("extdata/cpa_cov_", sub(" ", "", k), ".tif"), 
                                       format="GTiff", options="COMPRESS=LZW", overwrite=T); file.remove(cpa_files); rm(cpa_cov); gc()
    }
  }
})

dat <- raster::raster("extdata/cpa_cov_I-II.tif")
dat
rm(list=ls()); gc()

# Get coverage of mpa & cpa data
lapply(iucn_cat2, function(k){
  if(!file.exists(paste0("extdata/pa_cov_", sub(" ", "", k), ".tif"))){
    if(file.exists(paste0("extdata/WDPA_Dec2020_All_", sub(" ", "", k), ".rds"))){
      pa <- readRDS(paste0("extdata/WDPA_Dec2020_All_", sub(" ", "", k), ".rds"))
      
      # Specify raster resolution (1km)
      r <- raster::raster(nrow=21600, ncol=43200, xmn=-180, xmx=180, ymn=-90, ymx=90,
                          crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      r_sub <- raster::crop(r, as(pa, "Spatial"))
      tmpdir <- file.path(tempdir(), "splitRaster")
      dir.create(tmpdir)
      r_sub <- SpaDES.tools::splitRaster(r_sub, nx=5, ny=5, path=file.path(tmpdir, "mpa_sub"))
      lapply(1:length(r_sub), function(l){
        if(!file.exists(paste0("data/pa_cov_", sub(" ", "", k), "_", l, ".tif"))){
          pa_cov <- fasterize::fasterize(pa, r_sub[[l]])
          pa_cov <- raster::aggregate(pa_cov, fact=10, fun=sum, 
                                       filename=paste0("extdata/pa_cov_", sub(" ", "", k), "_", l, ".tif"), 
                                       format="GTiff", options="COMPRESS=LZW", overwrite=T) 
        }
      }); rm(pa, r_sub); gc()
    }
  }
})

lapply(iucn_cat2, function(k){
  if(!file.exists(paste0("extdata/pa_cov_", sub(" ", "", k), ".tif"))){
    pa_files <- list.files(path="extdata", pattern=paste0("pa_cov_", sub(" ", "", k), "_"), full.names=T)
    if(length(pa_files) == 25){
      pa_cov <- lapply(pa_files, raster)
      pa_cov <- do.call(raster::merge, c(pa_cov, tolerance=0.5))
      r <- raster::raster(nrow=2160, ncol=4320, xmn=-180, xmx=180, ymn=-90, ymx=90,
                          crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
      pa_cov <- raster::extend(pa_cov, r)
      pa_cov <- raster::projectRaster(pa_cov, crs="+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs", 
                                       filename=paste0("extdata/pa_cov_", sub(" ", "", k), ".tif"), 
                                       format="GTiff", options="COMPRESS=LZW", overwrite=T); file.remove(pa_files); rm(pa_cov); gc()
    }
  }
})

dat <- raster::raster("extdata/pa_cov_I-II.tif")
dat
rm(list=ls()); gc()

########################################