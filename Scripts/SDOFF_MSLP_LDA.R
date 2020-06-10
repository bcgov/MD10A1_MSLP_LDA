library(readr)
library(tidyverse)
library(candisc)
library(MASS)

#Read in SDoff cluster results from SDOFF_CLUSTER.R output
sdoff_clust<-read_csv("RawData/k_means_clusters.csv")

#Add year column
sdoff_clust$year<-c(2000:2018)

#Read Fall ERA5 PCA scores from FALL_ERA5_PCA.R output
era5_scores<-read_csv("RawData/ERA5_PCA_SCORES.csv")
era5_scores$X1<-NULL

lda_matx <- sdoff_clust %>% 
  dplyr::select(k_clust,year) %>%
  dplyr::left_join(era5_scores,by = 'year')

lda_matx<-as.data.frame(lda_matx)

lda_matx[,c(3:43)]<-scale(lda_matx[,c(3:43)])

init_pop<-function(N,max_terms,max_pc)
{
  pop<-list()
  
  for(e in c(1:N))
  {
    
    n_terms<-sample(c(1:max_terms),1)
    
    pc_vec<-c()
    
    for(i in c(1:n_terms))
    {
      pc_vec[i]<-sample(c(1:max_pc),1)
    }
    
    
    pop[[e]]<-pc_vec
  }
  
  return(pop)
}


vec_to_f<-function(pc_vec)
{
  f<-paste0("PC",pc_vec,collapse="+")
  
  f<-paste("lda_matx$k_clust ~",f)
  
  return(f)
}


calc_fit<-function(pc_vec,data,obs_clust)
{
  
  f<-vec_to_f(pc_vec)
  test <- lda(as.formula(f), data = data, CV=TRUE)
  
  mu<-mean(obs_clust == test$class)
  
  n<-length(obs_clust)
  p<-length(pc_vec)
  
  mu_adj<- 1-((n-1)/(n-p))*(1-mu)
  
  return(mu_adj)
}


init_fit<-function(pop)
{
  
  fit<-c()
  
  for(e in c(1:length(pop)))
  {
    fit[e]<-calc_fit(pop[[e]],lda_matx,lda_matx$k_clust)
  }
  
  return(fit)
}


select_parents<-function(prctl,fit,pop)
{
  thresh<-quantile(fit,.25)
  
  idx_vec<-c(1:length(fit))
  
  canidates<-idx_vec[fit>=thresh]
  
  P1<-pop[[sample(canidates,1)]]
  P2<-pop[[sample(canidates,1)]]
  
  return(list(P1,P2))
  
}


cross_over<-function(parents)
{
  P1<-parents[[1]]
  P2<-parents[[2]]
  comb<-unique(c(P1,P2))
  
  n<-length(comb)
  
  child<-sample(comb,sample(c(1:n),1))
  
  return(child)
}



mutation<-function(element,max_pc,mute_rate)
{
  size<-length(element)
  
  for(i in c(1:size))
  {
    rand<-runif(1)
    
    if(rand<mute_rate)
    {
      element[i]<-sample(c(1:max_pc),1)
    }
    
  }
  return(element)
}




evolve<-function(N,max_terms,max_pc,genz,prctl,mute_rate,data,obs_clust)
{
  
  pop<-init_pop(N,max_terms,max_pc)
  fit<-init_fit(pop)
  
  log_vec<-c(mean(fit))
  
  for(g in c(1:genz))
  {
    parents<-select_parents(prctl,fit,pop)
    
    child<-cross_over(parents)
    
    child<-mutation(child,max_pc,mute_rate)
    
    child_fit<-calc_fit(child,data,obs_clust)
    
    min_idx<-which.min(fit)
    min_fit<-fit[min_idx]
    
    if(child_fit>min_fit)
    {
      pop[[min_idx]]<-child
      fit[min_idx]<-child_fit
    }
    
    log_vec[g]<-mean(fit)
  }
  
  
  return(list(pop,fit,log_vec))
  
}

evolved<-evolve(1000,15,23,150000,.8,.15,lda_matx,lda_matx$k_clust)

plot(evolved[[3]])

best_mod<-evolved[[1]][[which.max(evolved[[2]])]]


lda_mod<-lda(f=as.formula(vec_to_f(best_mod)),data=lda_matx)

plot(lda_mod)

hm<-lapply(evolved[[1]],function(y) {sum(y)})



x=lm(cbind(PC1,PC2,PC6,PC13,PC22)~k_clust, lda_matx)
lda_can<-candisc(x, term="k_clust")
lda_scores<-lda_can$scores
lda_can$coeffs.std
