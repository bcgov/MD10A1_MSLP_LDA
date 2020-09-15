library(rsoi)
library(tidyverse)
library(tibble)
library(ratmos)

oni <- download_oni() %>%
  filter(Month=='Sep' | Month=='Oct' | Month=='Nov') %>%
  group_by(Year) %>%
  summarise(ONI = mean(ONI))

soi <- download_soi() %>%
  filter(Month=='Sep' | Month=='Oct' | Month=='Nov') %>%
  group_by(Year) %>%
  summarise(SOI = mean(SOI))

ao <- download_ao() %>%
  filter(Month=='Sep' | Month=='Oct' | Month=='Nov') %>%
  group_by(Year) %>%
  summarise(AO = mean(AO))


aao <- download_aao() %>%
  filter(Month=='Sep' | Month=='Oct' | Month=='Nov') %>%
  group_by(Year) %>%
  summarise(AAO = mean(AAO))


nao <- download_nao() %>%
  filter(Month=='Sep' | Month=='Oct' | Month=='Nov') %>%
  group_by(Year) %>%
  summarise(NAO= mean(NAO))

pna<- ratmos::get_pna() %>%
  filter(Month == 9 | Month == 10 | Month == 11) %>%
  group_by(Year) %>%
  summarise(PNA = mean(PNA))


pdo <- read_csv('RawData/pdo_tele.csv') %>%
  as_tibble() %>%
  mutate(PDO = (SEP+OCT+NOV)/3) %>%
  dplyr::rename(Year = YEAR) %>%
  select(Year,PDO)


tele <- as_tibble(c(1979:2018)) %>%
  dplyr::rename(Year = value)


tele <- tele %>%
  full_join(oni, by = 'Year') %>%
  full_join(soi, by = 'Year') %>%
  full_join(ao, by = 'Year') %>%
  full_join(aao, by = 'Year') %>%
  full_join(nao, by = 'Year') %>%
  full_join(pna, by = 'Year') %>%
  full_join(pdo, by = 'Year') %>%
  dplyr::filter(Year>=1979 & Year<=2018)


ERA5_PCA_SCORES <- read_csv("RawData/ERA5_PCA_SCORES.csv") %>%
  as_tibble() %>%
  dplyr::rename(Year = year) %>%
  dplyr::select(-one_of('X1'))


corr_tabl <- tele %>%
  full_join(ERA5_PCA_SCORES,by='Year')


write_csv(corr_tabl,"RawData/ERA5_TELE_CORRS.csv")

cor(corr_tabl, method = 'spearman', use = "complete.obs")


#Tests based on subset LDA results

cor.test(corr_tabl$ONI,corr_tabl$PC1, method = 'spearman')

cor.test(corr_tabl$SOI,corr_tabl$PC1, method = 'spearman', exact = T)

cor.test(corr_tabl$AAO,corr_tabl$PC2, method = 'spearman')

cor.test(corr_tabl$AO,corr_tabl$PC3, method = 'spearman')

cor.test(corr_tabl$PNA,corr_tabl$PC11, method = 'spearman')




