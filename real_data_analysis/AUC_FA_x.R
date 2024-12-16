library(tidyverse)
library(refund)
library(survival)
library(pROC)
library(mgcv)
library(orthogonalsplinebasis)

source('source_code/functions.R')
source('R_code/new/est_long.R')

load('RData/new/FPCs_new1.RData')
load('RData/new/cov_est_new.RData')

load('dataset/new/multi-omics_imputed.RData')
V = data_all$V
N = data_all$N
Y_all = data_all$Y
metadata_all = data_all$metadata
P_all = c(ncol(Y_all[[1]]), ncol(Y_all[[2]]))  ## 749, 104

## 
table(table(metadata_all[[1]]$ID))

F0 = 1
Fl = rep(NA, V)
v = 1
Fl[v] = 2  ## 40.3% total variation

v = 2
Fl[v] = 6  # 44.5% total variation


L0 = result_all$C0_result$L0
L = vector('list', V)
L[[1]] = c(result_all$cov_est_result[[1]]$L, L0)
L[[2]] = c(result_all$cov_est_result[[2]]$L, L0)


## Extract essential information
K = rep(NA, V)
time_all = vector('list', V)
ID_all = vector('list', V)
for (v in 1:V) {
  K[v] = nrow(Y_all[[v]])
  time = metadata_all[[v]]$time
  ID = metadata_all[[v]]$ID
  
  time_all[[v]] = time
  ID_all[[v]] = ID
}

scores = zeta_est

##
long = read.csv('dataset/processed/long.csv')
surv = read.csv('dataset/surv.csv')
surv = surv[surv$RID %in% long$RID, ]

long.bl = long[long$Years_bl==0, ]
all.equal(long.bl$RID, surv$RID)  ## TRUE

## Fit Cox model
c.test <- function(fit1, fit2){
  ctest = concordance(fit1, fit2)
  contr <- c(-1, 1)
  dtest <- contr %*% coef(ctest)
  dvar <- contr %*% vcov(ctest) %*% contr
  c(contrast=dtest, sd=sqrt(dvar), z=dtest/sqrt(dvar), p = 2*(1-pnorm(abs(dtest/sqrt(dvar))))) 
  p = 2*(1-pnorm(abs(dtest/sqrt(dvar))))
  return(p)
}

## FPC scores from lipidomics + gut_metabolite
surv2 = cbind(surv, long.bl$age, long.bl$sex, long.bl$APOE4, scores, xi_est[[1]], xi_est[[2]])
colnames(surv2)[5:ncol(surv2)] = c('Age', 'Sex', 'APOE4',  
                                   paste0('zeta_', 1:L0), 
                                   paste0('xi1_',  1:sum(L[[1]][1:Fl[1]])), 
                                   paste0('xi2_', 1:sum(L[[2]][1:Fl[2]])))
surv2 = surv2 %>% select(-RID, -ad_bl)
cox.obj = coxph(Surv(surv_time, status) ~ ., data = surv2)
concordance(cox.obj)  # 0.6971


## NULL model
cox.obj2 = coxph(Surv(surv_time, status) ~ Age + Sex + APOE4, data = surv2)
concordance(cox.obj2)  # 0.6771
c.test(cox.obj, cox.obj2)  ## 0.027

## Lipid only 
load('RData/new/FPCs_lipidomics_new.RData')

L1 = ncol(xi_est_lipid)
surv3 = cbind(surv, long.bl$age, long.bl$sex, long.bl$APOE4, xi_est_lipid)
colnames(surv3)[5:ncol(surv3)] = c('Age', 'Sex', 'APOE4', 
                                   paste0('xi_lipid_', 1:L1))
surv3 = surv3 %>% select(-RID, -ad_bl)
cox.obj.lipid = coxph(Surv(surv_time, status) ~ ., data = surv3)
concordance(cox.obj.lipid)  # 0.6901
c.test(cox.obj, cox.obj.lipid)  # 0.409

## Gut metabolite only 
load('RData/real_data_analysis/FPCs_gut_new.RData')

L1 = ncol(xi_est_gut)
surv3 = cbind(surv, long.bl$age, long.bl$sex, long.bl$APOE4, xi_est_gut)
colnames(surv3)[5:ncol(surv3)] = c('Age', 'Sex', 'APOE4',
                                   paste0('xi_gut_', 1:L1))
surv3 = surv3 %>% select(-RID, -ad_bl)
cox.obj.gut = coxph(Surv(surv_time, status) ~ ., data = surv3)
concordance(cox.obj.gut)  # 0.6912
c.test(cox.obj, cox.obj.gut)  # 0.415

## Baseline lipidomics + gut metabolomics ##
load("RData/real_data_analysis/FPCs_baseline_new.RData")
surv3 = cbind(surv, long.bl$age, long.bl$sex, long.bl$APOE4, Z_est[[1]], Z_est[[2]], U)
colnames(surv3)[5:ncol(surv3)] = c('Age', 'Sex', 'APOE4', paste0('Z1_', 1:Fl[1]), 
                                   paste0('Z2_', 1:Fl[2]), paste0('U_', 1:F0))
surv3 = surv3 %>% select(-RID, -ad_bl)
cox.obj.base = coxph(Surv(surv_time, status) ~ ., data = surv3)
concordance(cox.obj.base)  # 0.6906
c.test(cox.obj, cox.obj.base) # 0.312

library(xtable)
#save.image(file = 'RData/real_data_analysis/LS-MOFA.RData')

#load('RData/real_data_analysis/LS-MOFA.RData')
summ = summary(cox.obj)$coefficients
ci = summary(cox.obj)$conf.int

summ1 = cbind(summ[, c(2, 3, 5)], ci[, c(3, 4)])
summ1 = cbind(summ1, NA)
summ1[, c(1:3)] = round(summ1[, c(1:3)], 3)

summ1[, 6] = paste0(round(summ1[, 4], 3), ', ', round(summ1[, 5], 3))
summ1 = summ1[, c(-4, -5)]
xtable(summ1[, c(1, 2, 4, 3)])

