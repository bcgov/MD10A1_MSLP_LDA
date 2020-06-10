library(tidyverse)
library(bcdata)
library(bcmaps)
library(RColorBrewer)
library(factoextra)
library(cluster)

####Define Functions####

#Function for converting time series raster data into a table suitable for S-Mode PCA. *Assumes specific formating of tile names and folder structure. 
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

#Function returns a mean image of SDoff data for years specified by vector 'yr_vect'. SDoff band is 2, set -9 to NA. !!Assumes specific formatting of M*D10A1 tiff names, i.e., see global var snow_dir_path"!!
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
      
    }
    else
    {
      #Get the complete dir path as a string
      path<-paste(snow_dir_path,year,'.tif',sep = "")
      
      #Read in the snow duration raster for the given year
      rast <-raster(path, band=2)
      
      
      #Add the new year to the raster stack 
      s <- stack(s,rast)
    }
    i<-i+1
  }
  
  NAvalue(s)<--9.0
  
  s_mean<-calc(s,mean, na.rm=TRUE)
  
  return(s_mean)
}

#Function for wrting S-Mode PCA loadings to GeoTiff for import to GIS platforms...also returns raster object. Provide Lat Lon coordinates, a raster profile, prcomp (PCA) object, and the PC number (Eigenvalue index). 
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

#####Declair global varibles####

#Start and end over which M*D10A1 snow duration data is available
mod_start_year<-2000
mod_end_year<-2018

#Path to gridded M*10A1 derived snow duration Geotiff's (ommiting year)
snow_dir_path <-"RawData/Annual_Snow_Metrics/MD10A1_SD_"


#Get Eco_Prov as geometry, exclude 'NEP'
ecoprov <- ecoprovinces() %>% 
  filter(ECOPROVINCE_CODE!='NEP')

ecoprov<-ecoprov$geometry

#Get boundry of British Columbia as SF object
bc_boun<-bc_bound()$geometry

#Establish a color ramp, for sure the best one. 
cols <- brewer.pal(11, "Spectral")

##### Run Sdoff PCA and K-Means Cluster Analysis, generate mean Cluster plots. ####

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
plot(prop_varex, type = "b", xlab = "Principle Component",ylab = "Proportion of Variance Explained (%)")
dev.off()

#Get PCA scores
mod_scores <- as.data.frame(dur_pca$x)

#Use silhouette plot to compare kmeans with 2,3,4,5 or 6 clusters. Use PC scores as input to speed up computation. 
set.seed(101)

k_clust<-kmeans(x=mod_scores, centers=2,nstart = 1000)
sil<-silhouette(k_clust$cluster,dist(mod_scores))
fviz_silhouette(sil)

k_clust<-kmeans(x=mod_scores, centers=3,nstart = 1000)
sil<-silhouette(k_clust$cluster,dist(mod_scores))
fviz_silhouette(sil)

k_clust<-kmeans(x=mod_scores, centers=4,nstart = 1000)
sil<-silhouette(k_clust$cluster,dist(mod_scores))
fviz_silhouette(sil)

k_clust<-kmeans(x=mod_scores, centers=5,nstart = 1000)
sil<-silhouette(k_clust$cluster,dist(mod_scores))
fviz_silhouette(sil)

k_clust<-kmeans(x=mod_scores, centers=6,nstart = 1000)
sil<-silhouette(k_clust$cluster,dist(mod_scores))
fviz_silhouette(sil)

#!!!!Decided from silhouette plot qualitativley!!!!
K<-2

#Perform k-means cluster analysis on the snow duration PCA results using k chosen from above silhouette charts 
clusters <- kmeans(x=mod_scores, centers = K, nstart = 1000)


class_tabl<-as.data.frame(cbind(clusters$cluster,mod_scores))

names(class_tabl)[1]<-'k_clust'

write.csv(class_tabl,"RawData/k_means_clusters.csv")



