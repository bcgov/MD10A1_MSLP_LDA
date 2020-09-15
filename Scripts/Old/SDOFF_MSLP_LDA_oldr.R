library(readr)
library(tidyverse)
library(candisc)
library(MASS)

#Read in SDoff cluster results from SDOFF_CLUSTER.R output
sdoff_clust<-read_csv("RawData/sdoff_kcluster_results.csv")

#Add year column
sdoff_clust$year<-c(2000:2018)
sdoff_clust$X1<-NULL

#Read Fall ERA5 PCA scores from FALL_ERA5_PCA.R output
era5_scores<-read_csv("RawData/ERA5_PCA_SCORES.csv")
era5_scores$X1<-NULL

lda_matx <- sdoff_clust %>% 
  dplyr::select(k_clust,year) %>%
  dplyr::left_join(era5_scores,by = 'year')

lda_matx<-as.data.frame(lda_matx)


x=lm(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15)~k_clust, lda_matx)

lda_can<-candisc(x, term="k_clust", ndim=1)
lda_scores<-lda_can$scores

acc_vec<-c()

for(row in c(1:nrow(lda_matx)))
{
  x_loo=lm(cbind(PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15)~k_clust, lda_matx[-row,])
  
  lda_can_loo<-candisc(x_loo, term="k_clust", ndim=1)
  lda_scores_loo<-lda_can_loo$scores
  
  train_x<-lda_scores_loo[,2]
  train_y<-lda_scores_loo[,1]
  
  test_x<-lda_scores[row,2]
  test_y<-lda_scores[row,1]
  
  yhat<-knn(train = train_x,  test = test_x, cl=train_y,k=4,prob = TRUE)
  
  acc_vec[row]<-yhat
}



print(c("Leave-One-Out Acc.",mean(acc_vec)))

table(acc_vec,lda_matx$k_clust)


lda_can$coeffs.std

ggplot(lda_scores, aes(group=k_clust, y=Can1))+geom_boxplot()

Anova(x, test="Wilks")