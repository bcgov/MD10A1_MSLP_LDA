####Set global varibles####
library(readr)
library(tidyverse)
library(candisc)
library(MASS)
library(hrbrthemes)
library(xtable)

#Initial population size for GA
N<-2000

#Maximum number of terms for a LDA model
max_terms<-15

#Highest MSLP PC to concider for LDA model, based of 1% prop. varaince threshold 
max_pc<-15

#Number of GA generations 
generations<-1500000

#Percentile Adj. R2 threshold for fitness in GA selection 
prctl<-.75

#GA mutation rate 
mute_rate<-.15


####Generate Regression Table for LDA ####

#Read in SDoff cluster results from SDOFF_CLUSTER.R output
sdoff_clust<-read_csv("RawData/sdoff_kcluster_results.csv")

#Add year column
sdoff_clust$year<-c(2000:2018)

#Read Fall ERA5 PCA scores from FALL_ERA5_PCA.R output
era5_scores<-read_csv("RawData/ERA5_PCA_SCORES.csv")
era5_scores$X1<-NULL

#Join both tables by year
lda_matx <- sdoff_clust %>% 
  dplyr::select(k_clust,year) %>%
  dplyr::left_join(era5_scores,by = 'year')

#As a data frame 
lda_matx<-as.data.frame(lda_matx)


####Define Functions####

#Function randomly initlizes a subset of feature for LDA
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

#Converts a vector of feature index to a formula
vec_to_f<-function(pc_vec)
{
  f<-paste0("PC",pc_vec,collapse="+")
  
  f<-paste("lda_matx$k_clust ~",f)
  
  return(f)
}

#Calculates the fitness (Adj R2) of a LDA model provided as a PC vector (pc_vec) 
calc_fit<-function(pc_vec,data,obs_clust)
{
  
  f<-vec_to_f(pc_vec)
  test <- lda(as.formula(f), data = data, CV=TRUE)
  
  mu<-mean(obs_clust == test$class)
  
  #n<-length(obs_clust)
  #p<-length(pc_vec)
  
  #mu_adj<- 1-((n-1)/(n-p))*(1-mu)
  
  return(mu)
}

#Calculates the fitness of an intial population 
init_fit<-function(pop)
{
  
  fit<-c()
  
  for(e in c(1:length(pop)))
  {
    fit[e]<-calc_fit(pop[[e]],lda_matx,lda_matx$k_clust)
  }
  
  return(fit)
}

#Selects two parents randomly from members in the population with fitness above a specified percentile 
select_parents<-function(prctl,fit,pop)
{
  thresh<-quantile(fit,prctl)
  
  idx_vec<-c(1:length(fit))
  
  canidates<-idx_vec[fit>=thresh]
  
  P1<-pop[[sample(canidates,1)]]
  P2<-pop[[sample(canidates,1)]]
  
  return(list(P1,P2))
  
}

#Randomly combines two parents to create a offspring model 
cross_over<-function(parents)
{
  P1<-parents[[1]]
  P2<-parents[[2]]
  comb<-unique(c(P1,P2))
  
  n<-length(comb)
  
  child<-sample(comb,sample(c(1:n),1))
  
  return(child)
}

#Randomly mutate a element at a procided mutation rate (mute_rate)
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


#Function evolves a LDA model 
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

#Function for estimating the degree of seperation between two groups
calc_sep<-function(lda_mod,df)
{
  yhat<-predict(lda_mod,df)
  
  g1_mean<-mean(yhat$x[yhat$class==1])
  g1_var<-var(yhat$x[yhat$class==1])
  
  g2_mean<-mean(yhat$x[yhat$class==2])
  g2_var<-var(yhat$x[yhat$class==2])
  
  ratio<-(g1_mean-g2_mean)^2/(g1_var+g2_var)
  
  return(ratio)
  
}

#Function for plotting LD1 scores @ k=2
plot_multi_histogram <- function(df, feature, label_column) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black",bins = 19) +
    geom_density(alpha=0.7) +
    geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density") +
    theme_classic(base_size = 14)
  plt + guides(fill=guide_legend(title=expression(SD[OFF]~Cluster)))+theme(aspect.ratio = 1)
}

####Perform GA feature selection####

#Set seed for reproducability 
set.seed(101)

#Evolve a LDA model 
evolved<-evolve(N,max_terms,max_pc,generations,prctl,mute_rate,lda_matx,lda_matx$k_clust)

####Summarise LDA from GA Feature Selection####

#Get the vctor of Population average fitness and plot by generation 
mean_r2<-evolved[[3]]
mean_r2<-as.data.frame(cbind(mean_r2,c(1:length(mean_r2))))
colnames(mean_r2)[2]<-"Generation"

ggplot(mean_r2)+geom_line(aes(x=Generation,y=mean_r2)) + theme_classic() + labs(y = "Population Avg. Fitness")

ggsave('Manuscript/tatolatex/Figures/LDA/GA_PopAvgLOOCV_TS.jpeg')


sep_vec<-c()

for(i in c(1:length(evolved[[1]])))
{
  lda_temp_mod<-lda(f=as.formula(vec_to_f(evolved[[1]][[i]])),data=lda_matx)
  
  sep_vec[i]<-calc_sep(lda_temp_mod,lda_matx)
  
}

#Choose model that maximizes distance between group means, and minimizes the variance of each group, and has least number of inputs
mods<-as.data.frame(cbind(evolved[[1]],sep_vec))
mods<-mods[order(as.numeric(mods$sep_vec)),]

#Looks at result, choose best model
view(mods)

#Type in best model index below, (320)
best_mod<-evolved[[1]][[320]]

#Fit the LDA for the best model using LOO Cross-Validation
lda_mod<-lda(f=as.formula(vec_to_f(best_mod)),data=lda_matx,CV=T)

#Print the overall accuracy 
mean(lda_mod$class==lda_matx$k_clust)

#Print the LOOCV confusion matrix 
loocv<-table(lda_mod$class,lda_matx$k_clust)

caption<-"Leave-one-out cross-validate accuracy matrix for selected LDA model."

xtable(loocv,caption=caption,label = "tab:loocv")

anov.mod <- lm(cbind(PC1,PC2,PC3,PC8,PC10,PC11) ~ k_clust, data=lda_matx)
Anova(anov.mod, test="Wilks")


#Refit LDA model without LOOCV
lda_mod<-lda(f=as.formula(vec_to_f(best_mod)),data=lda_matx)

#Get model summary 
lda_mod

#Get unstand. coeff. 
Load_Coef_NS<-as.data.frame(t(lda_mod$scaling))

colnames(Load_Coef_NS)<-c('PC10','PC11','PC02','PC08','PC03','PC01')

Load_Coef_NS<-Load_Coef_NS[,order(colnames(Load_Coef_NS))]

caption<-"Descrimainat function coefficents for selected LDA model."

xtable(Load_Coef_NS,label="tab:lda_coef",caption=caption, digits = 5)

lda_group_means<-as.data.frame(lda_mod$means)

colnames(lda_group_means)<-c('PC01','PC10','PC11','PC02','PC03','PC08')

caption<-"The mean Fall ERA5 PC scores by sdoff{} cluster for all PC slected for input into LDA"

lda_group_means<-lda_group_means[,order(colnames(lda_group_means))]


xtable(lda_group_means,label="tab:lda_group_means",caption<-caption)

ggplot(lda_matx) + geom_boxplot( aes(x = as.factor(k_clust), y = PC1)) + theme_classic(base_size = 28) + labs(y = expression(MSLP[FALL]~PC1),x = expression(SD[OFF]~Cluster))
ggsave('Manuscript/tatolatex/Figures/LDA/PC1_Box.png')
ggplot(lda_matx) + geom_boxplot( aes(x = as.factor(k_clust), y = PC2)) + theme_classic(base_size = 28) + labs(y = expression(MSLP[FALL]~PC2),x = expression(SD[OFF]~Cluster))
ggsave('Manuscript/tatolatex/Figures/LDA/PC2_Box.png')
ggplot(lda_matx) + geom_boxplot( aes(x = as.factor(k_clust), y = PC3)) + theme_classic(base_size = 28) + labs(y = expression(MSLP[FALL]~PC3),x = expression(SD[OFF]~Cluster))
ggsave('Manuscript/tatolatex/Figures/LDA/PC3_Box.png')
ggplot(lda_matx) + geom_boxplot( aes(x = as.factor(k_clust), y = PC8)) + theme_classic(base_size = 28) + labs(y = expression(MSLP[FALL]~PC8),x = expression(SD[OFF]~Cluster))
ggsave('Manuscript/tatolatex/Figures/LDA/PC8_Box.png')
ggplot(lda_matx) + geom_boxplot( aes(x = as.factor(k_clust), y = PC10)) + theme_classic(base_size = 28) + labs(y = expression(MSLP[FALL]~PC10),x = expression(SD[OFF]~Cluster))
ggsave('Manuscript/tatolatex/Figures/LDA/PC10_Box.png')
ggplot(lda_matx) + geom_boxplot( aes(x = as.factor(k_clust), y = PC11)) + theme_classic(base_size = 28) + labs(y = expression(MSLP[FALL]~PC11),x = expression(SD[OFF]~Cluster))
ggsave('Manuscript/tatolatex/Figures/LDA/PC11_Box.png')


#Get LD1 scores and plot by class (score plot for k=2)
LD1_Scores<-predict(lda_mod)$x
LD1_Class<-predict(lda_mod)$class
LD1<-as.data.frame(cbind(LD1_Class,LD1_Scores))
LD1$LD1_Class<-as.factor(LD1_Class)


plot_multi_histogram(LD1,'LD1','LD1_Class')

ggsave('Manuscript/tatolatex/Figures/LDA/LD1_ScorePlt.png')


#Scale MSLP PC data 
lda_matx_scaled<-as.data.frame(scale(lda_matx[,c(3:ncol(lda_matx))]))

#Refit LDA on scaled data 
lda_mod<-lda(f=as.formula(vec_to_f(best_mod)),data=lda_matx_scaled)


#Plot standardized coefficents (k=2)
Load_Coef<-as.data.frame(lda_mod$scaling)
Load_Coef$PC<-rownames(Load_Coef)

ggplot(Load_Coef, aes(PC, LD1)) + geom_col(width=.25) + theme_classic() + labs(x="") + geom_hline(yintercept = 0)
ggsave('Manuscript/tatolatex/Figures/LDA/LD1_LoadingPlt.jpeg')


save.image(file = "SDoff_LDA_env.Rdata")
