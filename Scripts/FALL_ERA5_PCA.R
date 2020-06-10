library(reticulate)
library(stars)
library(tidyverse)
library(parallel)
library(RColorBrewer)


#Initilize python via reticulate
use_python("/usr/local/bin/python3")

#Call the ERA5 python download script 
source_python("Scripts/DL_Fall_ERA5_MSLP.py")

#Path to ERA5 download GRIB file 
era_path<-'RawData/download.grib'

#Initlize clusters
cl<-makeCluster(20)

#Function gets fall MSLP indicies within GRIB file for given year (!!Assumes data are from DL_Fall_ERA5_MSLP.py!!)
yr_to_grib_idx<-function(yr)
{
  strt_idx <- (yr*3)-5936
  end_idx <-  (yr*3)-5934
  
  return(c(strt_idx:end_idx))
}

#Function computes mean fall MSLP for given year 
agg_fall_era5<-function(yr)
{
  
  yr_rast<-read_stars(era_path,RasterIO = list(bands = yr_to_grib_idx(yr)))
  
  yr_rast<-st_apply(yr_rast, 1:2, mean, CLUSTER = cl) 
  
  yr_rast<-as.data.frame(yr_rast)
  
  return(yr_rast$mean)
}

#Initilize a fall_slp stack with first year (1979) 
fall_slp<-read_stars(era_path,RasterIO = list(bands = yr_to_grib_idx(1979)))
#Compute mean
fall_slp<-st_apply(fall_slp, 1:2, mean, CLUSTER = cl) 
fall_slp<-as.data.frame(fall_slp)
lon_lat<-fall_slp[,c(1,2)]
lon_lat$x[lon_lat$x>=180]<-lon_lat$x[lon_lat$x>=180]-360.0

fall_slp<-fall_slp$mean


for(yr in c(1980:2019))
{
  fall_slp<-rbind(fall_slp,agg_fall_era5(yr))
}

stopCluster(cl)

fall_slp<-as.data.frame(fall_slp)
slp_pca<-prcomp(fall_slp,scale.=T, center = T)


slp_scores<-as.data.frame(slp_pca$x)

slp_scores$year<-c(1979:2019)

write.csv(slp_scores,'RawData/ERA5_PCA_SCORES.csv')

loading<-rasterFromXYZ(cbind(lon_lat,slp_pca$rotation[,1]))
plot(loading)
cols <- brewer.pal(11, "Spectral")

coastlines <- st_read("ne-coastlines-10m/ne_10m_coastline.shp")$geometry

plot(loading, main = "",xlab = "",ylab="", box=FALSE,col = cols,xlim=c(-180,180), ylim=c(-90,90), asp=2)
plot(coastlines, add=T)


# hm<-varimax(slp_pca$rotation)
# loading<-rasterFromXYZ(cbind(lon_lat,hm$loadings[,3]))
# plot(loading)
