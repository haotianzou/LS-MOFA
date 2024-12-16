library(tidyverse)

data = read.csv('dataset/lipidomics/ADNI_covar.csv')
data = data %>% arrange(RID, timepoints)
length(unique(data$RID))

data_bl = data[data$timepoints=='bl', ]
table(data_bl$ad)
ID = data_bl$RID

## Filter out converter and ad_bl!=2
surv = read.csv('dataset/surv.csv')
surv2 = surv[surv$ad_bl!=2, ]
ID = surv2$RID

data = data[data$RID %in% ID, ]
data_bl = data[data$timepoints=='bl', ]  


## 
data_all = read.csv('dataset/lipidomics/ADNIMERGE_25May2023.csv')
data2 = data_all %>% mutate(timepoints = VISCODE) %>% select(RID, timepoints, Years_bl) %>% arrange(RID, Years_bl)
data2 = data2[data2$RID %in% ID, ]

data3 = data %>% inner_join(data2, by = c('RID', 'timepoints'))

lipid = read.csv('dataset/lipidomics/ADNI_BakerLipidome.csv')

## log 10 transformation and standardization
for (i in 2:782){
  tmp = lipid[, i]
  tmp = log10(tmp)
  tmp = unname(scale(tmp))[, 1]
  lipid[, i] = tmp
}

data4 = data3 %>% inner_join(lipid, by = 'RID_visit')

data5 = data4 %>% select(RID, timepoints, Years_bl, Sph.d18.1.:Total.OxSpecies)

data5_total_lipid = data5[, c(1:3, 785:831)]


## 781 lipid species
## Delete TG.SIM lipid and total lipid, we have 749 lipid species
data5.qc = data5 %>% select(-(TG.48.0...SIM.:TG.58.9...SIM.), -(Total.Sph:Total.OxSpecies))

## Filter the survival data 
data2 = data_all %>% mutate(timepoints = VISCODE) %>% select(RID, timepoints, Years_bl, DX) %>% arrange(RID, Years_bl)
data2 = data2[data2$RID %in% ID, ]
data2 = data2[data2$DX!="", ]
data2$ad = ifelse(data2$DX=='CN', 0, ifelse(data2$DX=='MCI', 1, ifelse(data2$DX=='Dementia', 2, NA)))
DX_last = data_bl$DXGrp_last_tp
data2$DXGrp_last_tp = rep(DX_last, table(data2$RID))

data3 = data2

surv = data.frame(RID = ID, surv_time = NA, status = NA, ad_bl = data_bl$ad)

for (i in 1:length(ID)){
  tmp.long = data3[data3$RID==ID[i], ]
  tmp.dx.last = tmp.long$DXGrp_last_tp[1]
  
  if (surv$ad_bl[i]==2) {
    surv$status[i] = 1
    surv$surv_time[i] = 0
  } else {
    ad_index = which(tmp.long$ad==2)
    if (length(ad_index)==0){
      surv$status[i] = 0
      surv$surv_time[i] = max(tmp.long$Years_bl)
    } else {
      surv$status[i] = 1
      surv$surv_time[i] = tmp.long$Years_bl[min(ad_index)]
    }
  }
  
  ## For those diagnoses varied across time. i.e., MCI, MCI, AD, MCI, MCI, MCI
  if (tmp.dx.last!='AD'){
    ad_index = which(tmp.long$ad==2)
    
    if (length(ad_index)!=0){
      if (min(ad_index)>1){
        surv$status[i] = 0
        surv$surv_time[i] = tmp.long$Years_bl[min(ad_index) - 1]
      }
    }
  }
  
}

my_ID = surv[surv$status==1 & surv$ad_bl!=2, ]$RID
tmp_ID = data_bl[data_bl$converter==1, ]$RID
all.equal(my_ID, tmp_ID)


## Sub: dataset without AD==2 at baseline 
## Long: dataset with AD==2 at baseline

surv_sub = surv[surv$ad_bl!=2, ]

AD.ID = surv[surv$ad_bl==2, ]$RID

data2 = data_all %>% mutate(timepoints = VISCODE) %>% select(RID, timepoints, Years_bl) %>% arrange(RID, Years_bl)
data2 = data2[data2$RID %in% ID, ]

data3 = data %>% inner_join(data2, by = c('RID', 'timepoints'))

long = data3 %>% arrange(RID, Years_bl)
write.csv(long, 'dataset/long.csv', row.names = F)

RID_visit_long = paste0(long$RID, "_", long$timepoints)

RID_visit = data3$RID_visit
lipid_long = data5.qc[data3$RID_visit %in% RID_visit_long, ]

write.csv(lipid_long, 'dataset/lipid_long.csv', row.names = F)

## 
ID.x = unique(long$RID)
max.followup = rep(NA, length(ID.x))
for (i in 1:length(ID.x)){
  tmp.long = long[long$RID==ID.x[i], ]
  max.followup[i] = max(tmp.long$Years_bl)
}

mean(max.followup)  ## 2.79
sd(max.followup) ## 1.83

## Join the baseline information
long = read.csv('dataset/long.csv')
lipid = read.csv('dataset/lipid_long.csv')

long.bl = long[which(long$timepoints=='bl'), ] %>% mutate(fasting = ifelse(fast=='Yes', 1, 0), 
                                                          cohort = ifelse(COHORT=='ADNI-2/GO', 1, 0), 
                                                          omega_3 = .OMEGA3_med, 
                                                          statin = .HYPERLIPIDAEMIA_med)
long.bl = long.bl %>% select(RID, age, sex, BMI, hdl, chol, trig, fasting, cohort, omega_3, statin, APOE4)

lipid.bl = lipid[lipid$timepoints=='bl', ] %>% select(-timepoints, -Years_bl)

full.data = surv %>% left_join(long.bl, by = 'RID') %>% left_join(lipid.bl, by = 'RID')
write.csv(full.data, 'dataset/full_data_baseline.csv', row.names = F)

## Read in gut metabolite
gut = read.csv('dataset/gut_metabolite/ADMCGUTMETABOLITESLONG_12_13_21_22Jul2023.csv')
gut = gut %>% arrange(RID, VISCODE2)
gut_bl = gut[gut$VISCODE2=='bl', ]
gut = gut[gut$RID %in% gut_bl$RID, ]

UID1 = data.frame(RID = unique(long$RID))
UID2 = data.frame(RID = unique(gut$RID))

UID_all = UID1 %>% inner_join(UID2, by = 'RID')

gut2 = gut[gut$RID %in% UID_all$RID, ]
lipid2 = lipid[lipid$RID %in% UID_all$RID, ]
long2 = long[long$RID %in% UID_all$RID, ]

## Extract the visit time for gut metabolite
library(lubridate)

UID = unique(gut2$RID)
Years_bl = NULL
for (i in 1:length(UID)){
  tmp.long = gut2[gut2$RID==UID[i], ]
  bl.date = ymd(tmp.long$EXAMDATE[1])
  date_diff = as.numeric(ymd(tmp.long$EXAMDATE) - rep(bl.date, nrow(tmp.long)))
  Years_bl = c(Years_bl, date_diff/365.25)
}

gut3 = gut2 %>% mutate(timepoints = VISCODE2, Years_bl = Years_bl) %>% select(RID, timepoints, Years_bl, 
                                                                              L_HISTIDINE:BUDCA)
## Log10 transform and scale
for (i in 4:ncol(gut3)){
  tmp = gut3[, i]
  tmp = ifelse(tmp==0, 0.1, tmp)
  tmp = log10(tmp)
  tmp = unname(scale(tmp))[, 1]
  gut3[, i] = tmp
}

write.csv(gut3, file = 'dataset/new/gut_metabolite.csv', row.names = F)
write.csv(lipid2, file = 'dataset/new/lipidomics.csv', row.names = F)
write.csv(long2, file = 'dataset/new/long.csv', row.names = F)

## Extract metadata
V = 2
data1 = read.csv('dataset/new/lipidomics.csv')  ## 3325 samples of lipidomics data
data2 = read.csv('dataset/new/gut_metabolite.csv') ## 3373 samples of gut metabolite
long = read.csv('dataset/new/long.csv') ## 3325 samples of longitudinal data
surv = read.csv('dataset/surv.csv')

surv = surv[surv$RID %in% long$RID, ]
uID = unique(long$RID)
N = length(unique(long$RID))  
data1.x = data1[0, ]
data2.x = data2[0, ]
long.x = long[0, ]

for (i in 1:N){
  surv.time = surv$surv_time[i]
  status = surv$status[i]
  
  tmp.data1 = data1[data1$RID==uID[i], ]
  tmp.data2 = data2[data2$RID==uID[i], ]
  tmp.data3 = long[long$RID==uID[i], ]
  
  if (status==1){
    tmp.data1 = tmp.data1[which(tmp.data1$Years_bl<surv.time), ]
    tmp.data2 = tmp.data2[which(tmp.data2$Years_bl<surv.time), ]
    tmp.data3 = tmp.data3[which(tmp.data3$Years_bl<surv.time), ]
  } else {
    tmp.data1 = tmp.data1[which(tmp.data1$Years_bl<=surv.time), ]
    tmp.data2 = tmp.data2[which(tmp.data2$Years_bl<=surv.time), ]
    tmp.data3 = tmp.data3[which(tmp.data3$Years_bl<=surv.time), ]
  }
  
  data1.x = rbind(data1.x, tmp.data1)
  data2.x = rbind(data2.x, tmp.data2)
  long.x = rbind(long.x, tmp.data3)
}

data1 = data1.x  ## 2813
data2 = data2.x  ## 2940
long = long.x

max.time = max(max(data1$Years_bl), max(data2$Years_bl))
N = length(unique(long$RID))  

long.bl = long[long$Years_bl==0, ]
table(long.bl$converter)

## 967 subjects with CN/MCI at baseline
## 671 non-converter (69.4%), 296 converter (30.6%)

Y_all = vector('list', V)
Y_all[[1]] = as.matrix(data1[, 4:ncol(data1)])
Y_all[[2]] = as.matrix(data2[, 4:ncol(data2)])

metadata_all = vector('list', V)
metadata_all[[1]] = data.frame(RID = data1$RID, ID = rep(1:N, table(data1$RID)), time = data1$Years_bl/max.time)
metadata_all[[2]] = data.frame(RID = data2$RID, ID = rep(1:N, table(data2$RID)), time = data2$Years_bl/max.time)
long$time = long$Years_bl/max.time

data_all = list(Y = Y_all, metadata = metadata_all, long = long, max.time = max.time, V = V, N = N)
save(list = 'data_all', file = 'dataset/new/multi-omics.RData')

tmp1 = table(data1$RID)
tmp2 = table(data2$RID)
index = which(tmp1!=tmp2)
df = data.frame(tmp1[index], tmp2[index])
df$surv.time = surv$surv_time[which(surv$RID %in% df$Var1)]
df$status = surv$status[which(surv$RID %in% df$Var1)]


d1 = data1[data1$RID==4714, ]
d2 = data2[data2$RID==4714, ]
