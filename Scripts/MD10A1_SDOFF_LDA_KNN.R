#Import required packages 
library(raster)
library(ncdf4)
library(doParallel)
library(parallel)
library(MASS)
library(RColorBrewer)
library(candisc)
library(bcmaps)
library(bcmapsdata)
library(MVN)
library(caret)
library(tidyverse)
library(readr)
library(corrplot)
library(ggplot2)
library(gridExtra)
library(ratmos)

#########################################################################################################
########################
###################################  Define Functions 
########################
#########################################################################################################

##(1)
#Function for converting time series raster data into a table suitable for S-Mode PCA. *Assumes specific formating of tile names and folder structure. 
##
rast_to_table<-function(strt_yr, end_yr, dir_string, b_num)
{
  
  #Get first year of snow melt date (Band 2 in image file) as a raster 
  s <- raster(paste(dir_string,strt_yr,'.tif',sep = ""), band=b_num)
  
  
  #Get important image metadata for generating loading images
  mod_rows<-nrow(s)
  mod_cols<-ncol(s)
  s_extent <- extent(s)
  s_rez <- res(s)
  s_crs <- crs(s)
  r <- raster(ncols=mod_cols, nrows=mod_rows)
  extent(r) <- s_extent
  res(r) <- s_rez
  crs(r) <- s_crs
  
  
  #For each year from start year+1 to end year, add the year to image stack 's'
  for (year in (strt_yr+1):end_yr){
    
    
    #Get the complete dir path as a string
    path<-paste(dir_string,year,'.tif',sep = "")
    
    #Read in the snow duration raster for the given year
    rast <-raster(path, band=b_num)
    
    #Add the new year to the raster stack 
    s <- stack(s,rast)
  }
  
  #Assighn -9 as NA value (Known From Meta Data)
  NAvalue(s) <- -9
  
  #Convert the raster data to a data frame, retain spatial coordinates 
  rast_data <-raster::as.data.frame(s, xy=TRUE)
  
  #Drop missing vectors
  rast_data <- rast_data[complete.cases(rast_data),]
  
  #Transpose the snow duration data frame, (i.e., columns = pixels, rows =  year) 
  rast_data <- t(rast_data)
  
  #Get logical of column indexes that have inter-annual varince, (i.e., not missing)
  has_var<-c(apply(rast_data[c(3:nrow(rast_data)),], 2, var) > 0)
  
  #Filter out cols(pixels) with varaince == 0
  rast_data<-rast_data[,has_var]
  
  #Get a data frame of the remaining Lat (C1) Lon (C2) indicies 
  lat_lon = as.data.frame(t(rast_data[c(1:2),]))
  
  #Get data frame without spatial coordinates columns (C>=3)
  rast_data<-rast_data[c(3:nrow(rast_data)),]
  
  #Return raster data, lat lon coords, and raster properties
  rtrn_lst<-list(rast_data,lat_lon,r)
  return(rtrn_lst)
}


#Function returns a mean image of SDoff data for years specified by vector 'yr_vect'.
#SDoff band is 2, set -9 to NA
#!!Assumes specific formatting of M*D10A1 tiff names, i.e., see global var snow_dir_path"!!
get_r_stack<-function(yr_vect,snow_dir_path)
{
  #Counter
  i<-0
  
  #Null raster stack
  s<-NULL
  
  #For each year add the tiff to raster stack s
  for (year in yr_vect)
  {
    if(i==0)
    {
      s <- raster(paste(snow_dir_path,year,'.tif',sep = ""), band=2)
      
      NAvalue(s) <- -9
    }
    else
    {
      #Get the complete dir path as a string
      path<-paste(snow_dir_path,year,'.tif',sep = "")
      
      #Read in the snow duration raster for the given year
      rast <-raster(path, band=2)
      
      NAvalue(s) <- -9
      
      #Add the new year to the raster stack 
      s <- stack(s,rast)
    }
    i<-i+1
  }
  
  
  s_mean<-calc(s,mean, na.rm=TRUE)
  
  return(s_mean)
}


#Function for wrting S-Mode PCA loadings to GeoTiff for import to GIS platforms...also returns raster object.
#Provide Lat Lon coordinates, a raster profile, prcomp (PCA) object, and the PC number (Eigenvalue index). 
grid_spat_load<-function(lat_lon,rast, pca, pc_num)
{
  #Get loading scores from PCA object at pc_num
  load<-c(pca$rotation[,pc_num])
  
  #Add corresponding spatial coords
  p<-data.frame(lat_lon, name=load)
  coordinates(p)<-~x+y
  
  #Convert to a raster object
  r<-rasterize(p,rast,'name',fun=mean)
  
  #Write out as tiff
  writeRaster(r, filename=paste("Manuscript/tatolatex/Figures/KNN/","PC",pc_num,"_loading.tif", sep=""), format="GTiff", overwrite=TRUE)
  
  #return raster object
  return(r)
}


#Function for getting mean annual MSLP for specified time period, converts spatial data to S-Mode PCA compatible matrix (rows = years, cols = pixels).
#Specify at start and end year to aggregate.
#Specify as start and end date to aggregate, e.g., '-09-01' and '-11-01' (Sep-Nov)
#Provide MSLP data as a nc object
aggr_slp<-function(slp_first_yr, slp_last_yr,slp_start_date,slp_end_date,slp)
{
  
  slp_mean<-NULL
  
  #For each year in range specified 
  for (yr in c(slp_first_yr:slp_last_yr))
  {
    
    #Establish sub annual time period over which to aggragate 
    srt<-paste(yr,slp_start_date,sep="")
    end<-paste(yr,slp_end_date,sep="")
    
    #Get a vector of time stamps whithin subannual time period specified 
    sub_perd<-time>=srt & time<=end
    
    #Get SLP data corresponding to sub annual time period 
    year<-slp[,,sub_perd]
    
    #Get number of months within the sub-annual time period 
    months<-dim(year)[3]
    
    #Get SLP data as a vector 
    slp_vec_long <- as.vector(year)
    
    #Each column as month, row as a pixel (long form) 
    slp_mat <- matrix(slp_vec_long, nrow=nlon*nlat, ncol=months)
    
    #Calculate row means, tranpose and create data frame 
    if (yr==slp_first_yr)
    {
      slp_mean<-as.data.frame(t(apply(slp_mat,1,mean)))
    }
    else
    {
      slp_mean[yr-slp_first_yr+1,]<-t(apply(slp_mat,1,mean))
    }
  }
  
  return(slp_mean)
  
}




#########################################################################################################
########################
###################################  Declair global varibles
########################
#########################################################################################################

#Start and end over which M*D10A1 snow duration data is available
mod_start_year<-2000
mod_end_year<-2018

#Path to gridded M*10A1 derived snow duration Geotiff's (ommiting year)
snow_dir_path <-"RawData/Annual_Snow_Metrics/MD10A1_SD_"

#Number of Sdoff clusters (!!From elbow plot analysis!!)
K<-3

#Path to gridded monthly Reanalysis 2 sea-level pressure
month_slp_pth<-"RawData/mslp.mon.mean.nc"

#Declair start and end year over which to compute principle components for MSLP data
slp_first_year=1979
slp_last_year=2018
slp_strt_date='-09-01'
slp_end_date='-11-01'

#Get Eco_Prov as geometry, exclude 'NEP'
ecoprov <- ecoprovinces() %>% 
  filter(ECOPROVINCE_CODE!='NEP')

ecoprov<-ecoprov$geometry

#Get boundry of British Columbia as SF object
bc_boun<-bc_bound()$geometry

#Establish a color ramp, for sure the best one. 
cols <- brewer.pal(11, "Spectral")


#########################################################################################################
########################
###################################  Run Sdoff PCA and K-Means Cluster Analysis, generate mean Cluster plots.
########################
#########################################################################################################

#Get the SDoff data into S-Mode matrix form, SDoff is band number 2
rast_data<-rast_to_table(mod_start_year,mod_end_year,snow_dir_path,2)

#Run PCA decompostion, rast_data is a list as -> [data table, lat_lon fields, empty raster object]
dur_pca <- prcomp(rast_data[[1]], center = TRUE, scale. = TRUE)

#Summarise the PCA results
summary(dur_pca)

#Generate a Scree-plot
std_dev <- dur_pca$sdev
pr_var <- std_dev^2
prop_varex <- (pr_var/sum(pr_var))*100
png("Manuscript/tatolatex/Figures/KNN/sdoff_scree.png")
plot(prop_varex, type = "h", xlab = "Principle Component",ylab = "Proportion of Variance Explained (%)")
dev.off()

#Get PCA scores
mod_scores <- as.data.frame(dur_pca$x)

#Generate Elbow plot to help decide on number of k clsuters, !!Update global vars!!
twss<-c()
for( c in c(1:10))
{
  cluster<-kmeans(x=mod_scores, centers=c, nstart = 10000)
  
  tss<-cluster$tot.withinss
  twss[c]<-tss
}

png("Manuscript/tatolatex/Figures/KNN/elbow_plot.png")
plot(twss, ylab = 'Total Within Sum of Squares', xlab='k-clusters', type = 'b')
abline(v=3)
dev.off()

#Perform k-means cluster analysis on the snow duration PCA results, *Must specify number of clusters to use in global vars*, drop last eignevalue (0%)
set.seed(20)
clusters <- kmeans(x=mod_scores[,c(1:18)], centers = K, nstart = 10000)


#View score averages by cluster (Two methods, same output)
clusters
mod_scores %>% 
  dplyr::mutate(clus=clusters$cluster) %>%
  dplyr::group_by(clus) %>%
  dplyr::summarise(pc1 = mean(PC1),pc2=mean(PC2),pc3=mean(PC3),pc4=mean(PC4))

#Generate year vector over period of record 
years<-c(mod_start_year:mod_end_year)

#Get years for each cluster type
clu_1<-years[clusters$cluster==1]
clu_2<-years[clusters$cluster==2]
clu_3<-years[clusters$cluster==3]

#Plot diffrence from time series mean for each SDoff Cluster 
ts_mean<-get_r_stack(years,snow_dir_path)
ts_mean<-projectRaster(ts_mean,crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
ts_mean[ts_mean<0]<--Inf
clu_1_r<-get_r_stack(clu_1,snow_dir_path)
clu_1_r<-projectRaster(clu_1_r,crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
clu_1_r[clu_1_r<0]<--Inf
clu_2_r<-get_r_stack(clu_2,snow_dir_path)
clu_2_r<-projectRaster(clu_2_r,crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
clu_2_r[clu_2_r<0]<--Inf
clu_3_r<-get_r_stack(clu_3,snow_dir_path)
clu_3_r<-projectRaster(clu_3_r,crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
clu_3_r[clu_3_r<0]<--Inf

clu_1_dif<-clu_1_r-ts_mean
clu_2_dif<-clu_2_r-ts_mean
clu_3_dif<-clu_3_r-ts_mean

clu_1_dif<-reclassify(clu_1_dif, rcl=c(-9999,-40,-50,
                                       -40,-30,-40,
                                       -30,-10,-30,
                                       -10,0,-10,
                                       0,10,10,
                                       10,30,30,
                                       30,40,40,
                                       40,9999,50))

clu_2_dif<-reclassify(clu_2_dif, rcl=c(-9999,-40,-50,
                                       -40,-30,-40,
                                       -30,-10,-30,
                                       -10,0,-10,
                                       0,10,10,
                                       10,30,30,
                                       30,40,40,
                                       40,9999,50))

clu_3_dif<-reclassify(clu_3_dif, rcl=c(-9999,-40,-50,
                                       -40,-30,-40,
                                       -30,-10,-30,
                                       -10,0,-10,
                                       0,10,10,
                                       10,30,30,
                                       30,40,40,
                                       40,9999,50))




RStoolbox::ggR(ts_mean, geom_raster = T) + 
  geom_sf(data = bc_boun, fill = NA, col = 'black' ) + 
  geom_sf(data=ecoprov, fill = NA, col = 'black') +
  coord_sf(xlim = c(1917109/9,1917109)) +
  scale_fill_gradientn(colours=cols, na.value = "white") +
  labs(x="",y="",fill="T.S. Mean (DSS)")+
  theme_void()+
  theme(legend.position  = c(.85,.7))

ggsave(filename = "Manuscript/tatolatex/Figures/KNN/ts_mean.jpeg", device = 'jpeg')




RStoolbox::ggR(clu_1_dif, geom_raster = T, forceCat = T) + 
  geom_sf(data = bc_boun, fill = NA, col = 'black' ) + 
  geom_sf(data=ecoprov, fill = NA, col = 'black') +
  coord_sf(xlim = c(1917109/9,1917109)) +
  scale_fill_manual(values = brewer.pal(10, "Spectral"), na.value="white") +
  labs(x="",y="",fill="Clust. 1 - T.S. Mean (Days)")+
  theme_void()+
  theme(legend.position  = c(.85,.7))

ggsave(filename = "Manuscript/tatolatex/Figures/KNN/cluster1_mean.jpeg", device = 'jpeg')

RStoolbox::ggR(clu_2_dif, geom_raster = T, forceCat = T) + 
  geom_sf(data = bc_boun, fill = NA, col = 'black' ) + 
  geom_sf(data=ecoprov, fill = NA, col = 'black') +
  coord_sf(xlim = c(1917109/9,1917109)) +
  scale_fill_manual(values = brewer.pal(10, "Spectral")) +
  labs(x="",y="",fill="Clust. 2 - T.S. Mean (Days)")+
  theme_void()+
  theme(legend.position  = c(.85,.7))

ggsave(filename = "Manuscript/tatolatex/Figures/KNN/cluster2_mean.jpeg", device = 'jpeg')

RStoolbox::ggR(clu_3_dif, geom_raster = T, forceCat = T) + 
  geom_sf(data = bc_boun, fill = NA, col = 'black' ) + 
  geom_sf(data=ecoprov, fill = NA, col = 'black') +
  coord_sf(xlim = c(1917109/9,1917109)) +
  scale_fill_manual(values = brewer.pal(10, "Spectral")) +
  labs(x="",y="",fill="Clust. 3 - T.S. Mean (Days)")+
  theme_void()+
  theme(legend.position  = c(.85,.7))

ggsave(filename = "Manuscript/tatolatex/Figures/KNN/cluster3_mean.jpeg", device = 'jpeg')


#Start parrallel cluster for generating loading images (Optional)
cl<-makeCluster(4)
registerDoParallel(cl)

#Run in Parrallel (!memory-hog!), get a list of first 4 spatialy reconstructed principle components. 
mod_load_imgs<-foreach(pc=1:4, .packages = 'raster') %dopar% + grid_spat_load(rast_data[[2]],rast_data[[3]],dur_pca,pc)


#Get loading raster from list as raster, reproject to EPSG:3005
pc1_ld<-mod_load_imgs[[1]][[1]]
pc1_ld<-projectRaster(pc1_ld,crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
pc2_ld<-mod_load_imgs[[2]][[1]]
pc2_ld<-projectRaster(pc2_ld,crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
pc3_ld<-mod_load_imgs[[3]][[1]]
pc3_ld<-projectRaster(pc3_ld,crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
pc4_ld<-mod_load_imgs[[4]][[1]]
pc4_ld<-projectRaster(pc4_ld,crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')

#Shutdown cluster object
stopCluster(cl)

# If only plotting: 
# pc1_ld <- raster("Manuscript/tatolatex/Figures/KNN/PC1_loading.tif")
# pc1_ld <- projectRaster(pc1_ld, crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
# pc2_ld <- raster("Manuscript/tatolatex/Figures/KNN/PC2_loading.tif")
# pc2_ld <- projectRaster(pc2_ld, crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
# pc3_ld <- raster("Manuscript/tatolatex/Figures/KNN/PC3_loading.tif")
# pc3_ld <- projectRaster(pc3_ld, crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
# pc4_ld <- raster("Manuscript/tatolatex/Figures/KNN/PC4_loading.tif")
# pc4_ld <- projectRaster(pc4_ld, crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')

RStoolbox::ggR(pc1_ld*10000, geom_raster = T, maxpixels = 1e6) + 
  geom_sf(data = bc_boun, fill = NA, col = 'black' ) +
  geom_sf(data=ecoprov, fill = NA, col = 'black') +
  coord_sf(xlim = c(1917109/9,1917109)) +
  scale_fill_gradientn(colours=cols, na.value = "white",
                       limits = c(-5,5),
                       breaks = seq(-5,5,2.5)) +
  guides(fill = guide_colorbar(barwidth = 2, barheight = 15, 
                               frame.colour = "black", ticks.colour = "black", 
                               draw.ulim = T, draw.llim = T, )) +
  labs(x="",y="",fill="Loading Coef. (x10000)")+
  theme_void() +
  theme(legend.position  = c(.85,.7))

ggsave(filename = "Manuscript/tatolatex/Figures/KNN/sdoff_pc1.jpeg", device = 'jpeg')

RStoolbox::ggR(pc2_ld*10000, geom_raster = T, maxpixels = 1e6) + 
  geom_sf(data = bc_boun, fill = NA, col = 'black' ) +
  geom_sf(data=ecoprov, fill = NA, col = 'black') +
  coord_sf(xlim = c(1917109/9,1917109)) +
  scale_fill_gradientn(colours=cols, na.value = "white",
                       limits = c(-10,10),
                       breaks = seq(-10,10,5)) +
  guides(fill = guide_colorbar(barwidth = 2, barheight = 15, 
                               frame.colour = "black", ticks.colour = "black", 
                               draw.ulim = T, draw.llim = T, )) +
  labs(x="",y="",fill="Loading Coef. (x10000)")+
  theme_void() +
  theme(legend.position  = c(.85,.7))

ggsave(filename = "Manuscript/tatolatex/Figures/KNN/sdoff_pc2.jpeg", device = 'jpeg')

RStoolbox::ggR(pc3_ld*10000, geom_raster = T, maxpixels = 1e6) + 
  geom_sf(data = bc_boun, fill = NA, col = 'black' ) +
  geom_sf(data=ecoprov, fill = NA, col = 'black') +
  coord_sf(xlim = c(1917109/9,1917109)) +
  scale_fill_gradientn(colours=cols, na.value = "white",
                       limits = c(-11,11),
                       breaks = seq(-10,10,5)) +
  guides(fill = guide_colorbar(barwidth = 2, barheight = 15, 
                               frame.colour = "black", ticks.colour = "black", 
                               draw.ulim = T, draw.llim = T, )) +
  labs(x="",y="",fill="Loading Coef. (x10000)")+
  theme_void() +
  theme(legend.position  = c(.85,.7))

ggsave(filename = "Manuscript/tatolatex/Figures/KNN/sdoff_pc3.jpeg", device = 'jpeg')


RStoolbox::ggR(pc4_ld*10000, geom_raster = T, maxpixels = 1e6) + 
  geom_sf(data = bc_boun, fill = NA, col = 'black' ) +
  geom_sf(data=ecoprov, fill = NA, col = 'black') +
  coord_sf(xlim = c(1917109/9,1917109)) +
  scale_fill_gradientn(colours=cols, na.value = "white",
                       limits = c(-15,15),
                       breaks = seq(-15,15,7.5)) +
  guides(fill = guide_colorbar(barwidth = 2, barheight = 15, 
                               frame.colour = "black", ticks.colour = "black", 
                               draw.ulim = T, draw.llim = T, )) +
  labs(x="",y="",fill="Loading Coef. (x10000)")+
  theme_void() +
  theme(legend.position  = c(.85,.7))

ggsave(filename = "Manuscript/tatolatex/Figures/KNN/sdoff_pc4.jpeg", device = 'jpeg')


# #Save loading rasters as JPEGS in manuscript folder
# jpeg("Manuscript/tatolatex/Figures/KNN/sdoff_pc1.jpeg", quality = 100)
# plot(pc1_ld, main="PC1", xlab="Easting", ylab = "Northing",col=cols)
# plot(ecoprov, add=TRUE)
# plot(bc_boun, add=TRUE)
# dev.off()
# jpeg("Manuscript/tatolatex/Figures/KNN/sdoff_pc2.jpeg",quality = 100)
# plot(pc2_ld, main="PC2", xlab="Easting", ylab = "Northing",col=cols)
# plot(ecoprov, add=TRUE)
# plot(bc_boun, add=TRUE)
# dev.off()
# jpeg("Manuscript/tatolatex/Figures/KNN/sdoff_pc3.jpeg",quality = 100)
# plot(pc3_ld, main="PC3", xlab="Easting", ylab = "Northing",col=cols)
# plot(ecoprov, add=TRUE)
# plot(bc_boun, add=TRUE)
# dev.off()
# jpeg("Manuscript/tatolatex/Figures/KNN/sdoff_pc4.jpeg",quality = 100)
# plot(pc4_ld, main="PC4", xlab="Easting", ylab = "Northing",col=cols)
# plot(ecoprov, add=TRUE)
# plot(bc_boun, add=TRUE)
# dev.off()



#Add a year column to the SDoff PCA scores table, represents hydrologic year, *assumes no missing years 
mod_scores$Year<- (mod_start_year-1) + seq(dim(mod_scores)[1])


#########################################################################################################
########################
################################### Read in global mean monthly sea level pressure derived from NCEP Reanalysis 2 data, run PCA on time series and plot spatial loadings.
########################
#########################################################################################################


#Read in the NCAR Reanlysis Monthly Mean SLP data 
ncin <- nc_open(month_slp_pth)

#print attributes to console 
print(ncin)

#Get each dimension within the ncd file 
lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon)
lat <- ncvar_get(ncin,"lat")
nlat <- dim(lat)
time <- ncvar_get(ncin,"time")
tunits <- ncatt_get(ncin,"time","units")
nt <- dim(time)


#Get the SLP Data 
slp <- ncvar_get(ncin,"mslp")

#Convert time vector to somthing R freindly 
time<-as.Date(time/24, origin='1800-01-01')

#Get mean MSLP for fall (Sep - Nov) period 
slp_mean<-aggr_slp(slp_first_year,slp_last_year,slp_strt_date,slp_end_date,slp)

#Run PCA on each year of mean SLP over the sub-annual period 
slp_pca <- prcomp(slp_mean, center = TRUE, scale. = TRUE)

#Summarise the MSLP PCA
summary(slp_pca)

#Generate a Scree-plot
std_dev <- slp_pca$sdev
pr_var <- std_dev^2
prop_varex <- (pr_var/sum(pr_var))*100
png("Manuscript/tatolatex/Figures/KNN/mslp_scree.png")
plot(prop_varex, type = "h", xlab = "Principle Component",ylab = "Proportion of Variance Explained (%)")
abline(h=prop_varex[15])
dev.off()

#Get the PCA scores as data frame 
slp_scores<-as.data.frame(slp_pca$x)

#Create a hydrologic year field 
slp_scores$Year<- (slp_first_year-1) + seq(dim(slp_scores)[1])

#Get observations during and after 2000 (MODIS record start)
sub_slp_scores<-slp_scores[slp_scores$Year>=mod_start_year,]



#########################################################################################################
########################
################################### Create final table for LDA anlysis, perform iterative model selection, and carry out descriptive LDA analysis 
########################
#########################################################################################################



#Add the sdoff cluster vector to MSLP PC score table 
reg_tab<-as.data.frame(cbind(sub_slp_scores, clusters$cluster))

#Convert to factor 
reg_tab$y_clust<- as.factor(reg_tab$`clusters$cluster`)
reg_tab$`clusters$cluster`<-NULL

#Attrib. to a responce vect y
y<-reg_tab$y_clust

#Check for multivariate normality of first 15 PC's
mvn(reg_tab[,c(1:15)])



reg_tab %>% 
  dplyr::group_by(y_clust) %>%
  dplyr::summarise_all(mean)



#Get raw and standardized LDA coef. using candisc package, *!Must manually enter model params! 
x=lm(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15)~y_clust, reg_tab)
lda_can<-candisc(x, term="y_clust")
lda_can
coef_std<-as.data.frame(lda_can$coeffs.std)
ld_scores<-as.data.frame(lda_can$scores)
lda_can$coeffs.raw
lda_can$coeffs.std


ld_scores %>% 
  dplyr::group_by(y_clust) %>%
  dplyr::summarise_all(mean)




png("Manuscript/tatolatex/Figures/KNN/scoreplot.png")
plot(Can2~Can1, data = as.data.frame(ld_scores), xlab = "LD1 Score", ylab="LD2 Score", main="LDA Score Plot", type = 'n')
text(ld_scores$Can1,ld_scores$Can2, labels = ld_scores$y_clust, pos = 3)
dev.off()


#Save loadingplot as JPEG
png("Manuscript/tatolatex/Figures/KNN/loadingplot.png")
plot(Can2~Can1, data = as.data.frame(coef_std), xlab = "LD1 Std. Coefficient", ylab="LD2 Std. Coefficient", main="LDA Loading Plot")
text(coef_std$Can1,coef_std$Can2, labels = row.names(coef_std), pos = 3)
abline(v=0, lty=2)
abline(h=0, lty=2)
dev.off()


#Asess leave-one-out classification accuracy of KNN model with 2 nearest neighbors. Expects specifc format of 'ld_scores' 
acc_vec<-c()

for(row in c(1:nrow(ld_scores)))
{
  
  yhat<-knn3Train(ld_scores[c(1:nrow(ld_scores))!=row,c(2,3)],ld_scores[row,c(2,3)],k=2,cl=y[c(1:nrow(ld_scores))!=row], prob = TRUE)
  
  
  acc_vec[row] = yhat==y[row]
}

print(c("Leave-One-Out Acc.",mean(acc_vec)))

table(acc_vec,y)

#########################################################################################################
########################
################################### Plot MSLP loading rasters for principle compoenntes included in final model (*from results above)
########################
#########################################################################################################



#Get PCA loadings 
slp_loadings <- slp_pca$rotation

#Manually specify PC vector, up to 15 PCs in this case
m_vect<-c(1:15)

#Generate a empty loading list 
loading_lst<-list()
cnt<-1


#Generate a list of MSLP loading rasters for each model MSLP PC input
for (ld in m_vect)
{
  
  pc_load <- slp_loadings[,ld]
  pc_load_mat <- matrix(pc_load, nrow = nlat, byrow = TRUE)
  pc_load_ras <- raster(pc_load_mat)
  extent(pc_load_ras) <- c(-180.0, 180.0, -90.0, 90.0)
  res(pc_load_ras) <- c((360.0/nlon),(180.0/nlat) )
  crs(pc_load_ras) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
  
  loading_lst[[cnt]]<-pc_load_ras
  cnt<-cnt+1
  
}


#Make plots of loading rasters in 'loading_lst' ...

#download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip", destfile = 'coastlines.zip')

#unzip(zipfile = "coastlines.zip", exdir = 'ne-coastlines-10m')

coastlines <- st_read("ne-coastlines-10m/ne_10m_coastline.shp")$geometry

for(i in c(1:15))
{
  
  path = paste0("Manuscript/tatolatex/Figures/KNN/mslp_pc",i,".jpeg")
  title = paste0("MSLP PC",i)
  
  jpeg(path,quality = 100)
  plot(loading_lst[[i]], main = "",xlab = "",ylab="", box=FALSE,col = cols,xlim=c(-180,180), ylim=c(-90,90), asp=2)
  plot(coastlines, add=TRUE)
  dev.off()
  
  
}




#########################################################################################################
########################
################################### Determine which of the MSLP PCs are correlated with primary teleconnections for discussion 
########################
#########################################################################################################


ao <- get_ao()
oni <- get_oni()
pna <- get_pna()
pdo <- get_pdo()
tele <- merge(merge(merge(oni,pdo),pna), ao) %>% 
  pivot_longer(c("ONI","PDO","PNA","AO"), names_to = "index", values_to = "value")


# Fall Indices
fall_oni <- tele %>%
  filter(index == 'ONI') %>%
  filter(Month %in% c(9,10,11)) %>% 
  group_by(Year, index) %>% 
  summarise(oni_fall = mean(value, na.rm = T))

fall_pdo <- tele %>%
  filter(index == 'PDO') %>%
  filter(Month %in% c(9,10,11)) %>% 
  group_by(Year, index) %>% 
  summarise(pdo_fall= mean(value, na.rm = T))

fall_pna <- tele %>%
  filter(index == 'PNA') %>%
  filter(Month %in% c(9,10,11)) %>% 
  group_by(Year, index) %>% 
  summarise(pna_fall = mean(value, na.rm = T))

fall_ao <- tele %>%
  filter(index == 'AO') %>%
  filter(Month %in% c(9,10,11)) %>% 
  group_by(Year, index) %>% 
  summarise(ao_fall = mean(value, na.rm = T))


mslp_retain<-paste('PC',c(1:15), sep="")



mslp_1_15<-reg_tab[,colnames(reg_tab)==mslp_retain]

mslp_1_15$Year<-c(2000:2018)



cor_tab<- mslp_1_15 %>%
  left_join(fall_oni, by = 'Year') %>%
  left_join(fall_pdo, by = 'Year') %>%
  left_join(fall_pna, by = 'Year') %>%
  left_join(fall_ao, by = 'Year')



write_csv(cor_tab, path = 'RawData/tele_mslp.csv')
