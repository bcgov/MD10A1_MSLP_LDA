library(tidyverse)
library(bcdata)
library(bcmaps)
library(RColorBrewer)
library(factoextra)
library(cluster)
library(doParallel)


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
  
  #Set the NA value of s to -9.0
  NAvalue(s)<--9.0
  
  s_mean<-calc(s,mean, na.rm=TRUE)
  
  return(s_mean)
}

#Function for wrting S-Mode PCA loadings to GeoTiff for import to GIS platforms...also returns raster object. Provide Lat Lon coordinates, a raster profile, prcomp (PCA) object, and the PC number (Eigenvalue index). 
grid_spat_load<-function(lat_lon,rast, pca, pc_num)
{
  #Get loading scores from PCA object at pc_num
  load<-c(pca$rotation[,pc_num])
  
  #Get as a XYZ table 
  xyz<-cbind(lat_lon,load)
  
  #Convert to a raster object
  r<-rasterFromXYZ(xyz,res=res(rast),crs=crs(rast))
  
  #Write out as tiff
  writeRaster(r, filename=paste("Manuscript/tatolatex/Figures/SDoff/","PC",pc_num,"_loading.tif", sep=""), format="GTiff", overwrite=TRUE)
  
  #return raster object
  return(r)
}

#Function to get loading raster from list as raster, reproject to EPSG:3005, create loading plot
plot_sdoff_load<-function(mod_load_imgs,num_load)
{
  #For each PC up to number of loadings (num_load)
  for(pc in c(1:num_load))
  {
    #Get the loading raster at PC and project 
    pc_ld<-mod_load_imgs[[pc]][[1]]
    pc_ld<-projectRaster(pc_ld,crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
    
    #Plot loading raster using RStoolbox
    RStoolbox::ggR(pc_ld*10000, geom_raster = T, maxpixels = 1e6) + 
      geom_sf(data = bc_boun, fill = NA, col = 'black' ) +
      geom_sf(data=ecoprov, fill = NA, col = 'black') +
      coord_sf(xlim = c(1917109/9,1917109)) +
      scale_fill_gradientn(colours=cols, na.value = "white") +
      guides(fill = guide_colorbar(barwidth = 1, barheight = 8, 
                                   frame.colour = "black", ticks.colour = "black", 
                                   draw.ulim = T, draw.llim = T, )) +
      labs(x="",y="",fill="Loading Coef.(x1e4)")+
      theme_void() +
      theme(legend.position  = c(.85,.7))
    
    #Save loading raster as JPEG
    out_file<-paste("Manuscript/tatolatex/Figures/SDoff/sdoff_pc",pc,".jpeg",sep="")
    ggsave(filename = out_file, device = 'jpeg')
  }
}

#Function for ploting cluster means given the corresponidng years as a vector, path to snow duration rasters, time series mean as raster, and the corresponding cluster number 
plot_clust_mean<-function(years,snow_dir_path,ts_mean,clust_num)
{
  #Get the mean raster for all years present in the years vector
  clust_r<-get_r_stack(years,snow_dir_path)
  clust_r<-projectRaster(clust_r,crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')
  
  #Diffrence the cluster mean from the time series mean 
  clust_dif<-clust_r-ts_mean
  
  #Get summary statistics for the diffrence raster for classifcation table 
  clust_smry<-summary(clust_dif)
  
  #Reclassify the diffrence raster 
  clust_dif<-reclassify(clust_dif, rcl=c(clust_smry[1],-30,-30,
                                         -30,-20,-25,
                                         -20,-10,-15,
                                         -10,-6,-8,
                                         -6,-2,-4,
                                         -2,2,0,
                                         2,6,4,
                                         6,10,8,
                                         10,20,15,
                                         20,30,25,
                                         30,clust_smry[5],30))
  
  #Generate output and label strings 
  labl_strg<-paste("Cluster ",clust_num,sep="")
  out_file<-paste("Manuscript/tatolatex/Figures/SDoff/cluster",clust_num,"_mean.jpeg",sep="")
  
  #Plot cluster diffrence anomoly using RStoolbox
  RStoolbox::ggR(clust_dif, geom_raster = T, forceCat = T) + 
    geom_sf(data = bc_boun, fill = NA, col = 'black' ) + 
    geom_sf(data=ecoprov, fill = NA, col = 'black') +
    coord_sf(xlim = c(1917109/9,1917109)) +
    scale_fill_manual(values = cols, na.value="white",labels = c("<-30","-30--20","-20--10","-10--6","-6--2","-2-2","2-6","6-10","10-20","20-30",">30")) +
    labs(x="",y="",fill=bquote(paste(.(labl_strg)," - Avg. ",SD[OFF]," (Days)")))+
    theme_void()+
    theme(legend.position  = c(.85,.7),legend.text=element_text(size=10))
  
  #Save as JPEG
  ggsave(filename = out_file, device = 'jpeg')
  
  
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

#!!!!Number of k-clusters, decided from silhouette plot qualitativley!!!!
K<-2

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
png("Manuscript/tatolatex/Figures/SDoff/sdoff_scree.png")
plot(prop_varex, type = "b", xlab = "Principle Component",ylab = "Proportion of Variance Explained (%)")
abline(h=1)
dev.off()


#Start parrallel cluster for generating loading images (Optional)
cl<-makeCluster(4)
registerDoParallel(cl)

#Run in Parrallel (!memory-hog!), get a list of first 4 spatialy reconstructed principle components. 
mod_load_imgs<-foreach(pc=1:18, .packages = 'raster') %dopar% + grid_spat_load(rast_data[[2]],rast_data[[3]],dur_pca,pc)

plot_sdoff_load(mod_load_imgs,18)

stopCluster(cl)

#Get PCA scores
mod_scores <- as.data.frame(dur_pca$x)

#Use silhouette plot to compare kmeans with 2,3,4,5 or 6 clusters. Use PC scores as input to speed up computation. 
set.seed(101)

clust_pal<-RColorBrewer::brewer.pal(4,'Set2')

k_clust<-kmeans(x=mod_scores, centers=2,nstart = 1000)
sil<-silhouette(k_clust$cluster,dist(mod_scores))
fviz_silhouette(sil,label = F, print.summary = T,palette = clust_pal,title="k = 2", legend.title="Cluster",font.main=c(16),font.y=c(16),font.legend=c(16),font.tickslab=c(14),ggtheme=theme_classic())+ylim(c(-.01,.3))
ggsave("Manuscript/tatolatex/Figures/SDoff/sil_k_2.png")

k_clust<-kmeans(x=mod_scores, centers=3,nstart = 1000)
sil<-silhouette(k_clust$cluster,dist(mod_scores))
fviz_silhouette(sil,label = F, print.summary = T,palette = clust_pal,title="k = 3", legend.title="Cluster",font.main=c(16),font.y=c(16),font.legend=c(16),font.tickslab=c(14),ggtheme=theme_classic())+ylim(c(-.01,.3))
ggsave("Manuscript/tatolatex/Figures/SDoff/sil_k_3.png")

k_clust<-kmeans(x=mod_scores, centers=4,nstart = 1000)
sil<-silhouette(k_clust$cluster,dist(mod_scores))
fviz_silhouette(sil,label = F, print.summary = T,palette = clust_pal,title="k = 4", legend.title="Cluster",font.main=c(16),font.y=c(16),font.legend=c(16),font.tickslab=c(14),ggtheme=theme_classic())+ylim(c(-.01,.3))
ggsave("Manuscript/tatolatex/Figures/SDoff/sil_k_4.png")

k_clust<-kmeans(x=mod_scores, centers=5,nstart = 1000)
sil<-silhouette(k_clust$cluster,dist(mod_scores))
fviz_silhouette(sil)

k_clust<-kmeans(x=mod_scores, centers=6,nstart = 1000)
sil<-silhouette(k_clust$cluster,dist(mod_scores))
fviz_silhouette(sil)

#Perform k-means cluster analysis on the snow duration PCA results using k chosen from above silhouette charts 

set.seed(101)
clusters <- kmeans(x=mod_scores, centers = K, nstart = 1000)

k_ts<-as.data.frame(t(cbind(c(2000:2018),clusters$cluster)))
colnames(k_ts)<-c(2000:2018)
k_ts<-k_ts[2,]

caption<-"Shown is the sdoff{} time series with the assoacated cluster value for each hydrologic year over the period 2002--2018."

xtable(k_ts,caption=caption, type="latex",label = "sdoff_k_ts")

class_tabl<-as.data.frame(cbind(clusters$cluster,mod_scores))

names(class_tabl)[1]<-'k_clust'

write.csv(class_tabl,"RawData/sdoff_kcluster_results.csv")


#Generate year vector over period of record 
years<-c(mod_start_year:mod_end_year)

#Get years for each cluster type
clu_1<-years[clusters$cluster==1]
clu_2<-years[clusters$cluster==2]
ts_mean<-get_r_stack(years,snow_dir_path)
ts_mean<-projectRaster(ts_mean,crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')

RStoolbox::ggR(ts_mean, geom_raster = T) + 
  geom_sf(data = bc_boun, fill = NA, col = 'black' ) + 
  geom_sf(data=ecoprov, fill = NA, col = 'black') +
  coord_sf(xlim = c(1917109/9,1917109)) +
  scale_fill_gradientn(colours=cols, na.value = "white") +
  labs(x="",y="",fill=expression(Average~SD[OFF]~(DSS)))+
  theme_void()+
  theme(legend.position  = c(.85,.7),legend.text=element_text(size=10))

ggsave(filename = "Manuscript/tatolatex/Figures/SDoff/ts_mean.jpeg", device = 'jpeg')

clu_1_r<-get_r_stack(clu_1,snow_dir_path)
clu_2_r<-get_r_stack(clu_2,snow_dir_path)
clu_dif<-clu_1_r-clu_2_r
summary(clu_dif)

clu_dif<-projectRaster(clu_dif,crs = '+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs')


#Get summary statistics for the diffrence raster for classifcation table 
clust_smry<-summary(clu_dif)


#Reclassify the diffrence raster 
clu_dif<-reclassify(clu_dif, rcl=c(clust_smry[1],-30,-30,
                                   -30,-20,-25,
                                   -20,-10,-15,
                                   -10,-6,-8,
                                   -6,-2,-4,
                                   -2,2,0,
                                   2,6,4,
                                   6,10,8,
                                   10,20,15,
                                   20,30,25,
                                   30,clust_smry[5],30))


#Plot cluster diffrence anomoly using RStoolbox
RStoolbox::ggR(clu_dif, geom_raster = T, forceCat = T) + 
  geom_sf(data = bc_boun, fill = NA, col = 'black' ) + 
  geom_sf(data=ecoprov, fill = NA, col = 'black') +
  coord_sf(xlim = c(1917109/9,1917109)) +
  scale_fill_manual(values = cols, na.value="white",labels = c("<-30","-30--20","-20--10","-10--6","-6--2","-2-2","2-6","6-10","10-20","20-30",">30")) +
  labs(x="",y="",fill="Clust. 1 - Clust. 2")+
  theme_void()+
  theme(legend.position  = c(.85,.7),legend.text=element_text(size=10))

ggsave(filename = "Manuscript/tatolatex/Figures/SDoff/clust_diff.jpeg", device = 'jpeg')
writeRaster(clu_dif,"Manuscript/tatolatex/Figures/SDoff/clust_diff.tif")



plot_clust_mean(clu_1,snow_dir_path,ts_mean,1)
plot_clust_mean(clu_2,snow_dir_path,ts_mean,2)


save.image(file = "SDoff_Cluster_env.Rdata")

