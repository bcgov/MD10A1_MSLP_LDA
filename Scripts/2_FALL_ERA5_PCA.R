####Global Vars####
library(reticulate)
library(stars)
library(tidyverse)
library(parallel)
library(RColorBrewer)
library(maptools)
library(raster)

colors <- brewer.pal(11, "Spectral")

coastlines <- as(st_read("ne-coastlines-10m/ne_10m_coastline.shp"),Class = "Spatial")

start_yr<-1979
end_yr<-2018

####Download ERA5 data using DL_Fall_ERA5_MSLP.py####

#Initilize python via reticulate
#use_python("/usr/local/bin/python3")

#Call the ERA5 python download script, if not called already, may have to edit output dir path in script.
#source_python("Scripts/DL_Fall_ERA5_MSLP.py")

#Path to ERA5 downloaded GRIB file 
era_path<-'RawData/fall_era5_download.grib'

####Define Functions####

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
  #Read the ERA5 GRIB file at bands corresponding to the year and months 
  yr_rast<-read_stars(era_path,RasterIO = list(bands = yr_to_grib_idx(yr)))
  
  #Reduce 3-months to seasonal (fall) mean 
  yr_rast<-st_apply(yr_rast, 1:2, mean, CLUSTER = cl) 
  
  #Convert to data frame
  yr_rast<-as.data.frame(yr_rast)
  
  #Return as vector 
  return(yr_rast$mean)
}

####Convert ERA5 monthly data into average fall MSLP in table form 

#Initlize clusters
cl<-makeCluster(20)

#Initilize a fall_slp stack with first year (1979) 
fall_slp<-read_stars(era_path,RasterIO = list(bands = yr_to_grib_idx(start_yr)))

#Compute mean
fall_slp<-st_apply(fall_slp, 1:2, mean, CLUSTER = cl) 
fall_slp<-as.data.frame(fall_slp)

#Get lat lon, convert to -180-180, -90-90
lon_lat<-fall_slp[,c(1,2)]
lon_lat$x[lon_lat$x>=180]<-lon_lat$x[lon_lat$x>=180]-360.0

#Get as vector 
fall_slp<-fall_slp$mean

#Add the remaining years to fall_slp table 
for(yr in c(start_yr+1:end_yr))
{
  fall_slp<-rbind(fall_slp,agg_fall_era5(yr))
}

#Shutdown clusters 
stopCluster(cl)

#Convert fall_slp to data frame and perform PCA
fall_slp<-as.data.frame(fall_slp)

####Run PCA on ERA5 fall MSLP data####

#Scale and center data, compute PCA,
slp_pca<-prcomp(fall_slp,scale.=T, center = T)

summary(slp_pca)

#Plot the screeplot 
std_dev <- slp_pca$sdev
pr_var <- std_dev^2
prop_varex <- (pr_var/sum(pr_var))*100

ggplot(as.data.frame(cbind(c(1:40),prop_varex)), aes(x=V1, y=prop_varex)) + 
  geom_point(size=3) + 
  geom_segment(aes(x=V1, 
                   xend=V1, 
                   y=0, 
                   yend=prop_varex)) + 
  labs(y="Proportion of Variance Explained (%)",x="Principle Component")+
  geom_hline(yintercept=prop_varex[15])+
  theme_classic()

ggsave("Manuscript/tatolatex/Figures/MSLP/mslp_scree.png")

#Get scores as dataframe 
slp_scores<-as.data.frame(slp_pca$x)

#Add year col to slp_scores 
slp_scores$year<-c(start_yr:end_yr)

#Write the ERA5 PCA scores to a csv
write.csv(slp_scores,'RawData/ERA5_PCA_SCORES.csv')

save.image(file = "ERA5_MSLP_env.Rdata")

####Plot MSLP loading rasters####

#Alternative plotting method, not implimented 
#spplot(loading,sp.layout=list('sp.lines', coastlines, lwd=2,first=F), scales = list(draw = TRUE))

for(pc in c(1:15))
{
  loading<-rasterFromXYZ(cbind(lon_lat,slp_pca$rotation[,pc]))
  crs(loading) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  
  path = paste0("Manuscript/tatolatex/Figures/MSLP/mslp_pc",pc,".jpeg")
  
  jpeg(path,quality = 100)
  plot(loading, main = "",xlab = "",ylab="",col=colors ,box=FALSE,xlim=c(-180,180), ylim=c(-90,90), asp=2)
  plot(coastlines, add = T)
  dev.off()
}

#Testing effect of varimax on loadings, not implimented 
#vari<-varimax(slp_pca$rotation)
#loading<-rasterFromXYZ(cbind(lon_lat,vari$loadings[,1]))
#plot(loading)
