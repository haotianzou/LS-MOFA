library(tidyverse)
library(gam)
library(mgcv)
library(refund)
library(matrixcalc)
library(splines)
library(orthogonalsplinebasis)
library(psych)
library(LaplacesDemon)
library(Matrix)
library(GPArotation)
library(survival)
library(tdROC)
library(ipred)
library(MOFA2)


N = 100
V = 2
P_all = rep(200, V)
Fl = rep(2, V)
F0 = 1
L = list(L1 = rep(2, Fl[[1]] + F0), L2 = rep(2, Fl[[2]] + F0))
Lu = list(L0 = rep(2, F0))

## Define the refined grid
t.argvals = seq(0, 1, by = 0.01)
lt = length(t.argvals)
uID_all = 1:N

source('R_code/sim_S3/data_gen_2omics_v3.R')
source('source_code/functions.R')

seed = 1

data = data_gen(N, V, P_all, Fl, F0, L, Lu, seed = 2024*seed)
Y_all = data$Y
metadata_all = data$metadata
W = data$W
surv.dat = data$surv.dat

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


## 
fname = paste0('result/sim_S3/LS-MOFA/', N, '/', seed, '.RData')
load(fname)

fname = paste0('result/sim_S3/LS-MOFA/', N, '/', seed, '_x.RData')
load(fname)

## Compute the estimated Y using LS-MOFA
U_est = result_all2$U_est
Z_est = result_all2$Z_est

W_est = vector('list', V)
lambda_est = vector('list', V)
for (v in 1:V){
  W_est[[v]] = result_all$cov_est_result[[v]]$W_est
  lambda_est[[v]] = result_all$Lambda.est[[v]]
}

Y_est = vector('list', V)
for (v in 1:V){
  tmp.U = matrix(unlist(U_est[[v]]), ncol = 1)
  tmp.Y = tmp.U %*% t(lambda_est[[v]]) + Z_est[[v]] %*% t(W_est[[v]])
  Y_est[[v]] = tmp.Y
}

## Compute the estimated Y using MEFISTO
outfile <- paste0("result/sim_S3/fit/", seed, '.hdf5')
fit = load_model(outfile, load_interpol_Z = T)
W2 = fit@expectations$W
Z2 = fit@expectations$Z
time2 = fit@covariates

Z2_name = names(Z2)
Z2_subj = as.numeric(substr(Z2_name, 6, str_length(Z2_name)))

W_MEFISTO = vector('list', V)
for (v in 1:V){
  W_MEFISTO[[v]] = matrix(NA, nrow = P_all[v], ncol = 5)
  
  tmp.W = W2[[v]] %>% as.data.frame()
  tmp.feature = rownames(tmp.W)
  tmp.locate = str_locate_all(tmp.feature, "_")
  
  tmp.feature_num = rep(NA, P_all[v])
  for (p in 1:P_all[v]){
    tmp.feature_num[p] = substr(tmp.feature[p], tmp.locate[[p]][1, 1] + 1, tmp.locate[[p]][2, 1] - 1)
  }
  tmp.feature_num = as.numeric(tmp.feature_num)
  
  tmp.W$feature = tmp.feature_num
  tmp.W = tmp.W %>% arrange(feature)
  W_MEFISTO[[v]] = as.matrix(tmp.W[, -6])
}

Y1_MEFISTO = vector('list', N)
Y2_MEFISTO = vector('list', N)
time_MEFISTO = vector('list', N)
intercept1_MEFISTO = vector('list', N)
intercept2_MEFISTO = vector('list', N)
for (i in 1:N){
  tmp.W1 = W_MEFISTO[[1]]
  tmp.W2 = W_MEFISTO[[2]]
  tmp.intercept1 = fit@intercepts[[1]]
  tmp.intercept2 = fit@intercepts[[2]]
  
  index = which(Z2_subj==i)
  tmp.Z = Z2[[index]]
  
  intercept1_MEFISTO[[i]] = tmp.intercept1[[index]]
  intercept2_MEFISTO[[i]] = tmp.intercept2[[index]]

  Y1_MEFISTO[[i]] = tmp.Z %*% t(tmp.W1) 
  Y2_MEFISTO[[i]] = tmp.Z %*% t(tmp.W2)
  time_MEFISTO[[i]] = time2[[index]]
}

theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

## Plot the true trajectory
f = function(view, subj, feature){
  
  ds = vector('list', V)
  for (v in 1:V) {
    tmp.Y = Y_all[[v]]
    tmp.metadata = metadata_all[[v]]
    
    tmp.Y = tmp.Y[which(tmp.metadata$ID==subj), ]
    tmp.time = tmp.metadata$time[which(tmp.metadata$ID==subj)]
    
    ds[[v]] = data.frame(ID = rep(subj, length(tmp.time)), 
                         time = tmp.time, Y = tmp.Y[, feature])
  }
  
  ds2 = vector('list', V)
  for (v in 1:V) {
    tmp.Y = Y_est[[v]]
    tmp.metadata = metadata_all[[v]]
    
    tmp.Y = tmp.Y[which(tmp.metadata$ID==subj), ]
    tmp.time = tmp.metadata$time[which(tmp.metadata$ID==subj)]
    
    ds2[[v]] = data.frame(ID = rep(subj, length(tmp.time)), 
                          time = tmp.time, Y = tmp.Y[, feature])
  }
  
  ds3 = vector('list', V)
  for (v in 1:V){
    if (v==1) tmp.Y = Y1_MEFISTO else tmp.Y = Y2_MEFISTO
    if (v==1) tmp.intercept = intercept1_MEFISTO else tmp.intercept = intercept2_MEFISTO
    
    tmp.Y = tmp.Y[[subj]]
    tmp.intercept = tmp.intercept[[subj]]
    tmp.metadata = metadata_all[[v]]
    
    tmp.time = tmp.metadata$time[which(tmp.metadata$ID==subj)]
    index = which(tmp.time %in% time_MEFISTO[[subj]][1, ])
    
    ds3[[v]] = data.frame(ID = rep(subj, length(tmp.time)), 
                          time = tmp.time, Y = tmp.Y[index, feature])
  }
  
  ds_all = rbind(ds[[view]], ds2[[view]], ds3[[view]])
  ds_all$label = rep(c('True', 'LS-MOFA', 'MEFISTO'), each = nrow(ds[[view]]))
  ds_all$subset = paste0('View ', view, ' Subj ', subj, ' Feature ', feature)
  
  return(ds_all)
}

ds_all1 = f(1, 1, 2)
ds_all2 = f(1, 1, 25)
ds_all3 = f(1, 20, 28)
ds_all4 = f(1, 80, 148)

ds_all5 = f(2, 1, 17)
ds_all6 = f(2, 1, 187)
ds_all7 = f(2, 20, 22)
ds_all8 = f(2, 80, 14)

ds_all = rbind(ds_all1, ds_all2, ds_all3, ds_all4, 
               ds_all5, ds_all6, ds_all7, ds_all8)

cairo_ps(filename = 'result/sim.eps', width = 16, height = 8)
ggplot(data = ds_all, aes(x = time, y = Y)) + 
  geom_point() + 
  geom_line(aes(linetype = label)) + 
  scale_linetype_manual(values = c('dashed', 'dotted', 'solid')) + 
  theme(legend.title = element_blank()) + 
  ylab('Omics value') + 
  facet_wrap(~ subset, nrow = 2, ncol = 4, scales='free')
dev.off()

ggsave(filename = 'result/sim.png', width = 16, height = 8)

