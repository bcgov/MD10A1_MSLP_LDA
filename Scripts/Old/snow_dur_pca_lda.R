
#Import required packages 
library(raster)
library(ggplot2)
library(ncdf4)
library(doParallel)
library(nnet)
library(parallel)
library(MASS)
library(caret)
library(bestNormalize)


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

#Start parrallel cluster for generating loading images (Optional)
numCores <- detectCores()
cl<-makeCluster(4)
registerDoParallel(cl)



#Get first year of snow duration (Band 3 in image file) as a raster 
s <- raster(paste(snow_dir_path,mod_start_year,'.tif',sep = ""), band=2)


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
for (year in (mod_start_year+1):mod_end_year){
  
  print(year)
  
  #Get the complete dir path as a string
  path<-paste(snow_dir_path,year,'.tif',sep = "")
  
  #Read in the snow duration raster for the given year
  rast <-raster(path, band=2)
  
  #Add the new year to the raster stack 
  s <- stack(s,rast)
}

#Remove temporary varibles
remove(path, year, rast) 

#NAvalue(s) <- -9

#Convert the raster data to a data frame, retain spatial coordinates 
rast_data <-raster::as.data.frame(s, xy=TRUE)

#Remove raster stack to save on RAM 
remove(s)

#Get a logical of complete rows 
non_miss<-complete.cases(rast_data)

#Drop missing rows 
rast_data<-rast_data[non_miss, ]

#Save RAM 
remove(non_miss)

#Scale only the snow duration columns, skip lat-lon
rast_data<-cbind(rast_data[,c(1:2)],scale(rast_data[,c(3:21)]))


#Transpose the snow duration data frame, (i.e., columns = pixels, rows =  year) 
rast_data <- t(rast_data)

#Get logical of column indexes that have inter-annual varince, (i.e., not missing)
has_var<-c(apply(rast_data[c(3:nrow(rast_data)),], 2, var) > 0.0)


#Filter out cols(pixels) with varaince == 0
rast_data<-rast_data[,has_var]

#Save on RAM 
remove(has_var)


#Get a data frame of the remaining Lat Lon indicies 
lat_lon = as.data.frame(t(rast_data[c(1:2),]))

#Get data frame without spatial coordinates columns 
rast_data<-rast_data[c(3:nrow(rast_data)),]


#Run PCA decompostion 
dur_pca <- prcomp(rast_data, center = TRUE, scale. = TRUE)

#Save on RAM
remove(rast_data)

#Summarise the PCA results
summary(dur_pca)

plot(dur_pca)


#Get PCA scores
mod_scores <- as.data.frame(dur_pca$x)

#write.csv(mod_scores, "C:\\Users\\hgleason\\Dropbox\\mod_tab.csv")

#Hart test function to determine rough number of 'k' clusters to use 
hart_test<-function(x_)
{
  for (i in c(1:ncol(x_)))
  {
    k<-i
    
    set.seed(20)
    clust_1 <- kmeans(x=x_, centers=i, nstart = 5000)
    set.seed(20)
    clust_2 <- kmeans(x=x_, centers=i+1, nstart = 5000)
    
    l<-clust_1$tot.withinss/clust_2$tot.withinss
    r<-nrow(x_)-k-1
    tot<-l*r
    
    print(c(k,tot))
    
    
  }
}

#Run hart test on mod_scores
hart_test(mod_scores)

#Perform k-means cluster analysis on the snow duration PCA results 
set.seed(20)
clusters <- kmeans(x=mod_scores, centers = 3, nstart = 5000)


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
years<-c(2000:2018)


summary(mod_scores$PC2)

#Get years for each cluster type
clu_1<-years[clusters$cluster==1]
clu_2<-years[clusters$cluster==2]
clu_3<-years[clusters$cluster==3]
# clu_4<-years[clusters$cluster==4]
# clu_5<-years[clusters$cluster==5]
# clu_6<-years[clusters$cluster==6]

#pc1_low<-years[mod_scores$PC2<=-741]
#pc1_hig<-years[mod_scores$PC2>=790]

#Get mean image for each cluster
clu_1_r<-get_r_stack(clu_1,snow_dir_path,1)
clu_2_r<-get_r_stack(clu_2,snow_dir_path,1)
clu_3_r<-get_r_stack(clu_3,snow_dir_path,1)
# clu_4_r<-get_r_stack(clu_4,snow_dir_path,1)
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
plot(clu_1_r)
plot(clu_2_r)
plot(clu_3_r)
# plot(clu_4_r)
# plot(clu_5_r)
# plot(clu_6_r)


clu_comb<-stack(clu_1_r,clu_2_r,clu_3_r)



plot(clu_comb, main= c("Cluster 1 Mean","Cluster 2 Mean","Cluster 3 Mean"),col = grey(1:255/255), box=FALSE)



# #Function for wrting loading to CSV with x,y,z fields for import to QGIS
grid_spat_load<-function(lat_lon,rast, pca, pc_num)
{
  load<-c(pca$rotation[,pc_num])
  p<-data.frame(lat_lon, name=load)
  coordinates(p)<-~x+y
  r<-rasterize(p,rast,'name',fun=mean)
  #writeRaster(r, filename=paste("PC",pc_num,"_loading.tif", sep=""), format="GTiff", overwrite=TRUE)
  return(r)
}


mod_load_imgs<-foreach(pc=1:4, .packages = 'raster') %dopar% + grid_spat_load(lat_lon,r,dur_pca,pc)
  

pc1_ld<-mod_load_imgs[[1]][[1]]
pc2_ld<-mod_load_imgs[[2]][[1]]
pc3_ld<-mod_load_imgs[[3]][[1]]
pc4_ld<-mod_load_imgs[[4]][[1]]

pc_ld_stk<-stack(pc1_ld, pc2_ld,pc3_ld,pc4_ld)


plot(pc_ld_stk, main = c("PC1","PC2","PC3","PC4"), col = colorRampPalette(c("black", "grey", "blue"))(255) )

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

#Print dimensions 
print(ncin$dim)

#Get the SLP Data 
slp <- ncvar_get(ncin,"mslp")

#Convert time vector to somthing R freindly 
time<-as.Date(time/24, origin='1800-01-01')


seasons<-list(list("-03-01","-05-01"),list("-06-01","-08-01"),list("-09-01","-11-01"),list("-12-01","-02-01"))


step_model<-function(s,seas_lst, slp, rmode)
{
  #Declair start and end date over which to compute median SLP
  slp_strt_date=seasons[[s]][[1]]
  slp_end_date=seasons[[s]][[2]]
  
  
  aggr_slp<-function(slp_first_yr, slp_last_yr,slp_start_date,slp_end_date,slp)
  {
    
    slp_median<-NULL
    
    #For each year in range specified 
    for (yr in c(slp_first_yr:slp_last_yr))
    {
      
      #Establish sub annual time period over which to aggragate 
      srt<-paste(yr,slp_start_date,sep="")
      end<-paste(yr,slp_end_date,sep="")
      
      if(s>3)
      {
        end<-paste(yr+1,slp_end_date,sep="")
      }
      
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
        slp_median<-as.data.frame(t(apply(slp_mat,1,mean)))
      }
      
      slp_median[yr-slp_first_yr+1,]<-t(apply(slp_mat,1,mean))
      
    }
    
    return(slp_median)
    
  }
  
  slp_median<-aggr_slp(slp_first_year,slp_last_year,slp_strt_date,slp_end_date,slp)
  
  
  #Run PCA on each year of mean SLP over the sub-annual period 
  slp_pca <- prcomp(slp_median, center = TRUE, scale. = TRUE)
  
  if(rmode==2)
  {
    return(slp_pca)
  }
  
  
  #Get the PCA scores as data frame 
  slp_scores<-as.data.frame(slp_pca$x)
  
  
  trans_list<-list()
  
  for (pc_ix in c(1:(length(slp_scores[1,])-1)))
  {
    pc<- paste("PC",pc_ix, sep="")
    
    pc_data<-slp_scores[pc]
    
    pc_data = as.matrix(as.data.frame(lapply(pc_data, as.numeric)))
    
    best_norm<-bestNormalize(pc_data,loo=TRUE)
    
    slp_scores[pc]<-best_norm$x.t
    
    trans_list[[pc_ix]]<-best_norm$chosen_transform
    
  }
  
  
  if(rmode==3)
  {
    return(slp_scores)
  }
  
  #Create a year field 
  slp_scores$Year<- (slp_first_year-1) + seq(dim(slp_scores)[1])
  
  #Get observations during and after 2002 (MODIS record start)
  sub_slp_scores<-slp_scores[slp_scores$Year>=mod_start_year,]
  
  
  
  reg_tab<-cbind(sub_slp_scores, clusters$cluster)
  
  if(rmode==1)
  {
    return(reg_tab)
  }
  
  
  std_dev <- slp_pca$sdev
  pr_var <- std_dev^2
  prop_varex <- pr_var/sum(pr_var)
  
  
  retain<-length(prop_varex[prop_varex>=0.01])
  
  
  
  best_mod_scr<-0.0
  best_mod<-NULL
  best_idx<-NULL
  
  for (e in c(1:retain))
  {
    for (g in c(e:retain))
    {
      for(h in c(g:retain))
      {
        for(i in c(h:retain))
        {
          for (j in c(i:retain))
          {
            for (k in c(j:retain))
            {
              
              if( e!=g & e!=h & e!=i & e!=j &  e!=k & g!=h & g!=i & g!=j & g!=k & h!=i & h!=j & h!=k & i!=j & i!=k & j!=k)
              {
                
                f<-paste0("PC",e,"+","PC",g,"+","PC",h,"+","PC",i,"+","PC",j,"+","PC",k, collapse = "+")
                f<-paste("reg_tab$`clusters$cluster` ~",f)
                
                
                test <- lda(as.formula(f), data = reg_tab, CV=TRUE)
                
                
                mu<-mean(reg_tab$`clusters$cluster`==test$class)
                
                if(mu>best_mod_scr)
                {
                  best_mod_scr<-mu
                  best_mod<-f
                  
                }
              }
            }
          }
        }
      }
    }
  }
  
  
  model<-lda(as.formula(best_mod),data = reg_tab, CV=TRUE)
  
  error<-c(s,mean(reg_tab$`clusters$cluster` == model$class))
  
  rtrn_lst<-list(model,error,trans_list)
  
  return(rtrn_lst)
  
  
}



cl<-makeCluster(4)
registerDoParallel(cl)


#Fit a LDA for each season specified in the season list 
model_lst<-foreach(s=c(1:4), .packages = c('MASS','bestNormalize')) %dopar% step_model(s,seasons,slp,0)
stopCluster(cl)

#Print the LOOV results for each seaon 
for (i in c(1:4))
{
  sub_lst<-model_lst[[i]]
  print(sub_lst[[2]])
  
}


#Pick season for final LDA anlysis
seas_idx<-2

#Get the stepwise model result for the specified season 
best_mod<-model_lst[[seas_idx]][[1]]

#If already run step-wise 
#best_mod<-"reg_tab$`clusters$cluster` ~ PC3 + PC13 + PC15 + PC19 + PC22 + PC29"

#Get the regression tabel assocated with the chosen season
reg_tab<-step_model(seas_idx,seasons,slp,1)

#Fit and print the LDA model results 
model<-lda(as.formula(best_mod),data = reg_tab)
model
plot(model)

#Generate confusionMatrix
model<-lda(as.formula(best_mod),data = reg_tab, CV=TRUE)
ob<-reg_tab$`clusters$cluster`
yhat<-model$class
table(ob,yhat)


#Get the SLP PCA for the specified season 
slp_pca<-step_model(seas_idx,seas_lst, slp,2)


#Get the PCA scores as data frame 
slp_scores<-as.data.frame(slp_pca$x)

#Create a year field 
slp_scores$Year<- (slp_first_year-1) + seq(dim(slp_scores)[1])

#Get observations during and after 2002 (MODIS record start)
sub_slp_scores<-slp_scores[slp_scores$Year>=mod_start_year,]


#Summarise results of SLP PCA
summary(slp_pca)


#Generate a Scree-plot
std_dev <- slp_pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
plot(prop_varex, type = "b", xlab = "Principle Component",ylab = "Proportion of Variance Explained")

ggplot( data=NULL, aes(x=c(1:39), y=cumsum(prop_varex))) +
  geom_line(linetype="dashed")+
  geom_point(size = 3) + 
  theme_bw() +
  labs(x = "Principle Component")+
  labs(y = "Cumulative Proportion of Variance Explained")


#Get PCA loadings 
slp_loadings <- slp_pca$rotation

#Get a specific loading, reshape to orignal geospatial extent 
loading_lst<-list()
cnt<-1

for (ld in c(1,4,5,8,22,23))
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


library(sf)

#download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip", destfile = 'coastlines.zip')

#unzip(zipfile = "coastlines.zip", exdir = 'ne-coastlines-10m')

coastlines <- st_read("ne-coastlines-10m/ne_10m_coastline.shp")$geometry

par(mfrow=c(3,2)) 

plot(loading_lst[[1]], main = "PC1", box=FALSE,col = grey(1:255/255))
plot(coastlines, add=TRUE)
plot(loading_lst[[2]], main = "PC4",box=FALSE,col = grey(1:255/255))
plot(coastlines, add=TRUE)
plot(loading_lst[[3]], main = "PC5",box=FALSE,col = grey(1:255/255))
plot(coastlines, add=TRUE)
plot(loading_lst[[4]], main = "PC8",box=FALSE,col = grey(1:255/255))
plot(coastlines, add=TRUE)
plot(loading_lst[[5]], main = "PC22",box=FALSE,col = grey(1:255/255))
plot(coastlines, add=TRUE)
plot(loading_lst[[6]], main = "PC23",box=FALSE,col = grey(1:255/255))
plot(coastlines, add=TRUE)



#View Loading on Map
plot(pc_load_ras)
mapview::mapView(pc_load_ras)
out<-'C:\\Users\\hgleason\\Downloads\\PC5.tif'
writeRaster(loading_lst[[3]], filename=out, format="GTiff", overwrite=TRUE)




