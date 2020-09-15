
#Import required packages 
library(raster)
library(ggplot2)
library(ncdf4)
library(doParallel)
library(parallel)
library(MASS)
library(caret)
library(sf)
library(tidyverse)

######################################################################################################################################################
######################################################################################################################################################
#Read in gridded annual snow duration over BC derived from M*D10A1 data (see prior study), run PCA on time series and export spatial loadings as CSV##
######################################################################################################################################################
######################################################################################################################################################

#########################
#Declair global varibles#
#########################

#Start and end over which M*D10A1 snow duration data is available
mod_start_year<-2000
mod_end_year<-2018

#Path to gridded M*10A1 derived snow duration Geotiff's (ommiting year)
snow_dir_path <-"/home/huntergleason/Dropbox/FLNRO/Projects/HYDAT/Basin_Scale/MODIS/Derived/Annual_Snow_Metrics/MD10A1_SD_"

#Path to gridded monthly Reanalysis 2 sea-level pressure
month_slp_pth<-"/home/huntergleason/Dropbox/FLNRO/Projects/LDA_MOD_SnowDur/RawData/mslp.mon.mean.nc"

#Declair start and end year over which to compute principle components
slp_first_year=1979
slp_last_year=2018


rast_to_table<-function(strt_yr, end_yr, dir_string, b_num)
{
  
  #Get first year of snow duration (Band 3 in image file) as a raster 
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
  
  
  
  #For each year 2003-2017, add the year to image stack 's'
  for (year in (strt_yr+1):end_yr){
    
    
    #Get the complete dir path as a string
    path<-paste(dir_string,year,'.tif',sep = "")
    
    #Read in the snow duration raster for the given year
    rast <-raster(path, band=b_num)
    
    #Add the new year to the raster stack 
    s <- stack(s,rast)
  }
  
  
  #NAvalue(s) <- -9
  
  #Convert the raster data to a data frame, retain spatial coordinates 
  rast_data <-raster::as.data.frame(s, xy=TRUE)
  
  #rast_data[,c(3:ncol(rast_data))]<-rast_data[,c(3:ncol(rast_data))]/366
  
  
  #rast_data <- rast_data[complete.cases(rast_data),]
  
  
  #Transpose the snow duration data frame, (i.e., columns = pixels, rows =  year) 
  rast_data <- t(rast_data)
  
  
  #Get logical of column indexes that have inter-annual varince, (i.e., not missing)
  has_var<-c(apply(rast_data[c(3:nrow(rast_data)),], 2, var) > 0)
  
  
  #Filter out cols(pixels) with varaince == 0
  rast_data<-rast_data[,has_var]
  
  rast_data[rast_data == -9] <- 1
  
  
  #Get a data frame of the remaining Lat Lon indicies 
  lat_lon = as.data.frame(t(rast_data[c(1:2),]))
  
  #Get data frame without spatial coordinates columns 
  rast_data<-rast_data[c(3:nrow(rast_data)),]
  
  
  rtrn_lst<-list(rast_data,lat_lon,r)
  
  return(rtrn_lst)
  
}

rast_data<-rast_to_table(2000,2018,snow_dir_path,2)


#Run PCA decompostion 
dur_pca <- prcomp(rast_data[[1]], center = TRUE, scale. = TRUE)


#Summarise the PCA results
summary(dur_pca)

#Generate a Scree-plot
std_dev <- dur_pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
plot(prop_varex, type = "b", xlab = "Principle Component",ylab = "Proportion of Variance Explained")

# ggplot( data=NULL, aes(x=c(1:length(prop_varex)), y=cumsum(prop_varex))) +
#   geom_line(linetype="dashed")+
#   geom_point(size = 3) + 
#   theme_bw() +
#   labs(x = "Principle Component")+
#   labs(y = "Cumulative Proportion of Variance Explained")


#Get PCA scores
mod_scores <- as.data.frame(dur_pca$x)

#Perform k-means cluster analysis on the snow duration PCA results 
set.seed(20)
clusters <- kmeans(x=mod_scores[,c(1:18)], centers = 3, nstart = 10000)

clusters


#Returns a mean image of snowduration for years specified by yr_vect
get_r_stack<-function(yr_vect,snow_dir_path,mode)
{
  i<-0
  s<-NULL
  for (year in yr_vect)
  {
    if(i==0)
    {
      s <- raster(paste(snow_dir_path,year,'.tif',sep = ""), band=2)
      
      s[s == -9] <- 1
    }
    else
    {
      #Get the complete dir path as a string
      path<-paste(snow_dir_path,year,'.tif',sep = "")
      
      #Read in the snow duration raster for the given year
      rast <-raster(path, band=2)
      
      rast[rast == -9] <- 1
      
      #Add the new year to the raster stack 
      s <- stack(s,rast)
    }
  }
  
  if(mode==1)
  {
    s_mean<-calc(s,mean)
  }
  else
  {
    s_mean<-calc(s,sd)
  }
  
  
  return(s_mean)
}

#Generate year vector over period of record 
years<-c(mod_start_year:mod_end_year)


#Get years for each cluster type
clu_1<-years[clusters$cluster==1]
clu_2<-years[clusters$cluster==2]
clu_3<-years[clusters$cluster==3]
#clu_4<-years[clusters$cluster==4]
# clu_5<-years[clusters$cluster==5]
# clu_6<-years[clusters$cluster==6]

#pc1_low<-years[mod_scores$PC2<=-741]
#pc1_hig<-years[mod_scores$PC2>=790]

#Get mean image for each cluster
clu_1_r<-get_r_stack(clu_1,snow_dir_path,1)
clu_2_r<-get_r_stack(clu_2,snow_dir_path,1)
clu_3_r<-get_r_stack(clu_3,snow_dir_path,1)
#clu_4_r<-get_r_stack(clu_4,snow_dir_path,1)
# clu_5_r<-get_r_stack(clu_5,snow_dir_path,1)
# clu_6_r<-get_r_stack(clu_6,snow_dir_path,1)
#clu_1_sd<-get_r_stack(clu_1,snow_dir_path,2)
#clu_2_sd<-get_r_stack(clu_2,snow_dir_path,2)
#clu_3_sd<-get_r_stack(clu_3,snow_dir_path,2)

#pc1_low_r<-get_r_stack(pc1_low,snow_dir_path)
#pc1_hig_r<-get_r_stack(pc1_hig,snow_dir_path)

#plot(pc1_low_r)
#plot(pc1_hig_r)

#Plot image for each cluster type 
plot(clu_1_r, main="Clust 1")
plot(clu_2_r, main="Clust 2")
plot(clu_3_r, main="Clust 3")
#plot(clu_4_r)
# plot(clu_5_r)
# plot(clu_6_r)


clu_comb<-stack(clu_1_r,clu_2_r,clu_3_r)



plot(clu_comb, main= c("Cluster 1 Mean","Cluster 2 Mean","Cluster 3 Mean"), box=FALSE)



# #Function for wrting loading to CSV with x,y,z fields for import to QGIS
grid_spat_load<-function(lat_lon,rast, pca, pc_num)
{
  load<-c(pca$rotation[,pc_num])
  p<-data.frame(lat_lon, name=load)
  coordinates(p)<-~x+y
  r<-rasterize(p,rast,'name',fun=mean)
  writeRaster(r, filename=paste("PC",pc_num,"_loading.tif", sep=""), format="GTiff", overwrite=TRUE)
  return(r)
}


#Start parrallel cluster for generating loading images (Optional)
cl<-makeCluster(4)
registerDoParallel(cl)


mod_load_imgs<-foreach(pc=1:4, .packages = 'raster') %dopar% + grid_spat_load(rast_data[[2]],rast_data[[3]],dur_pca,pc)


pc1_ld<-mod_load_imgs[[1]][[1]]
pc2_ld<-mod_load_imgs[[2]][[1]]
pc3_ld<-mod_load_imgs[[3]][[1]]
pc4_ld<-mod_load_imgs[[4]][[1]]

pc_ld_stk<-stack(pc1_ld, pc2_ld,pc3_ld,pc4_ld)


plot(pc_ld_stk, main = c("PC1","PC2","PC3","PC4"),frame = FALSE)


stopCluster(cl)

#Add a year column 
mod_scores$Year<- (mod_start_year-1) + seq(dim(mod_scores)[1])




##########################################################################################################################################################
##########################################################################################################################################################
#Read in gridded mean monthly sea level pressure derived from NCEP Reanalysis 2 data (see prior study), run PCA on time series and plot spatial loadings.#
##########################################################################################################################################################
##########################################################################################################################################################



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
    
    #Get number of month within the sub-annual time period 
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
    
    slp_mean[yr-slp_first_yr+1,]<-t(apply(slp_mat,1,mean))
    
  }
  
  return(slp_mean)
  
}

slp_mean<-aggr_slp(slp_first_year,slp_last_year,"-06-01","-08-01",slp)


#Run PCA on each year of mean SLP over the sub-annual period 
slp_pca <- prcomp(slp_mean, center = TRUE, scale. = TRUE)


#Get the PCA scores as data frame 
slp_scores<-as.data.frame(slp_pca$x)


#Create a year field 
slp_scores$Year<- (slp_first_year-1) + seq(dim(slp_scores)[1])

#Get observations during and after 2002 (MODIS record start)
sub_slp_scores<-slp_scores[slp_scores$Year>=mod_start_year,]


reg_tab<-as.data.frame(cbind(sub_slp_scores, clusters$cluster))

reg_tab$`clusters$cluster`<- as.factor(reg_tab$`clusters$cluster`)


retain<-length(prop_varex[prop_varex>=0.005])


best_mod_scr<-0.0
best_mod<-NULL
best_idx<-NULL


retain<-summary(slp_pca)$importance[2,]>=0.01

y_clust<-reg_tab$`clusters$cluster`

inputs<-reg_tab[,c(1:26)]
inputs<-cbind(inputs,y_clust)

library(candisc)

x=lm(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21,PC22,PC23,PC24,PC25,PC26)~y_clust,data=inputs)
out2=candisc(x, term='y_clust')

for( i in c(2:length(inputs)))
{
  
  test <- lda(x=inputs[,c(2,4,5,6,8,9,11,13,16)],grouping=y_clust, CV=TRUE)
  
  mu<-mean(y_clust == test$class)
  
  print(c(i,mu))
  
}



best_mod<-step_lda(matrix(reg_tab[,c(1:26)]),y_clust)


model<-lda(as.formula(best_mod),data = reg_tab, CV=TRUE)

error<-mean(reg_tab$`clusters$cluster` == model$class)

error

#Fit and print the LDA model results 
model<-lda(as.formula(best_mod),data = reg_tab)
model
plot(model, main = "LDA Score Plot")

#Generate confusionMatrix
model<-lda(as.formula(best_mod),data = reg_tab, CV=TRUE)
ob<-reg_tab$`clusters$cluster`
yhat<-model$class
table(ob,yhat)


#Summarise results of SLP PCA
summary(slp_pca)


#Generate a Scree-plot
std_dev <- slp_pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
plot(prop_varex, type = "b", xlab = "Principle Component",ylab = "Proportion of Variance Explained")



#Get PCA loadings 
slp_loadings <- slp_pca$rotation

#Get a specific loading, reshape to orignal geospatial extent 
loading_lst<-list()
cnt<-1

for (ld in c(1,5,11,12,14,17))
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




#download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip", destfile = 'coastlines.zip')

#unzip(zipfile = "coastlines.zip", exdir = 'ne-coastlines-10m')

coastlines <- st_read("ne-coastlines-10m/ne_10m_coastline.shp")$geometry

par(mfrow=c(3,2)) 

plot(loading_lst[[1]], main = "PC1", box=FALSE,col = grey(1:255/255))
plot(coastlines, add=TRUE)
plot(loading_lst[[2]], main = "PC5",box=FALSE,col = grey(1:255/255))
plot(coastlines, add=TRUE)
plot(loading_lst[[3]], main = "PC11",box=FALSE,col = grey(1:255/255))
plot(coastlines, add=TRUE)
plot(loading_lst[[4]], main = "PC12",box=FALSE,col = grey(1:255/255))
plot(coastlines, add=TRUE)
plot(loading_lst[[5]], main = "PC14",box=FALSE,col = grey(1:255/255))
plot(coastlines, add=TRUE)
plot(loading_lst[[6]], main = "PC17",box=FALSE,col = grey(1:255/255))
plot(coastlines, add=TRUE)



#View Loading on Map
# plot(pc_load_ras)
# mapview::mapView(pc_load_ras)
out<-'/home/huntergleason/Dropbox/FLNRO/Projects/LDA_MOD_SnowDur/Figures/slp_loadings/SLP_PC17.tif'
writeRaster(loading_lst[[6]], filename=out, format="GTiff", overwrite=TRUE)


mean_pc5pc12<-mean(loading_lst[[2]],loading_lst[[4]])
mean_pc11pc17<-mean(loading_lst[[3]],loading_lst[[6]])

diff<-mean_pc5pc12-mean_pc11pc17

plot(abs(diff),box=FALSE, main = "LD1 Sea Level Pressure Contrast")
plot(coastlines, add=TRUE)
pol <- rasterToPolygons(abs(diff), fun=function(x){x>quantile(abs(diff),.95)}, dissolve = TRUE)
plot(pol, add=TRUE)

diff<-loading_lst[[1]]-loading_lst[[5]]

plot(abs(diff),box=FALSE, main = "LD2 Sea Level Pressure Contrast")
plot(coastlines, add=TRUE)
pol <- rasterToPolygons(abs(diff), fun=function(x){x>quantile(abs(diff),.95)}, dissolve = TRUE)
plot(pol, add=TRUE)
