library(tidyverse)
library(psych)
library(LaplacesDemon)
library(GPArotation)
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
n.repeats = 100

MSE_all = vector('list', n.repeats)

## For each seed, compute MSE of Y and Y_hat

for (seed in 1:100){


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

##  Compute MSE of Y_est and Y_all
MSE1 = vector('list', V)
for (v in 1:V){
  tmp.MSE = (Y_est[[v]] - Y_all[[v]])^2
  MSE1[[v]] = mean(tmp.MSE)
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

Y_MEFISTO = vector('list', V)
for (v in 1:V){
  Y_MEFISTO[[v]] = matrix(NA, nrow = 0, ncol = P_all[v])
  time_v = time_all[[v]]
  metadata_v = metadata_all[[v]]
  if (v==1) tmp.Y = Y1_MEFISTO else tmp.Y = Y2_MEFISTO
  
  for (i in 1:N){
    tmp.time = as.vector(time_MEFISTO[[i]])
    tmp.Yx = tmp.Y[[i]]
    time_i = time_v[which(metadata_v$ID==i)]
    tmp.Y_sel = matrix(tmp.Yx[which(time_i %in% tmp.time), ], ncol = P_all[v])
    Y_MEFISTO[[v]] = rbind(Y_MEFISTO[[v]], tmp.Y_sel)
  }
}

##  Compute MSE of Y_est and Y_all
MSE2 = vector('list', V)
for (v in 1:V){
  tmp.MSE = (Y_MEFISTO[[v]] - Y_all[[v]])^2
  MSE2[[v]] = mean(tmp.MSE)
}

MSE_all[[seed]] = list(MSE1 = MSE1, MSE2 = MSE2)
}

save(list = c('MSE_all'), file = 'result/sim_S3_MSE.RData')

## Compare MSE
MSE_LSMOFA = MSE_MEFISTO = rep(NA, V)
SE_LSMOFA = SE_MEFISTO = rep(NA, V)
p = rep(NA, V)
for (v in 1:V){
  tmp.MSE1 = tmp.MSE2 = rep(NA, n.repeats)
  for (i in 1:n.repeats){
    tmp.MSE1[i] = MSE_all[[i]]$MSE1[[v]]
    tmp.MSE2[i] = MSE_all[[i]]$MSE2[[v]]
  }
  MSE_LSMOFA[v] = mean(tmp.MSE1)
  SE_LSMOFA[v] = sd(tmp.MSE1)
  
  MSE_MEFISTO[v] = mean(tmp.MSE2)
  SE_MEFISTO[v] = sd(tmp.MSE2)
  p[v] = wilcox.test(tmp.MSE1, tmp.MSE2, paired = T)$p.value
}

save(list = c('MSE_LSMOFA', 'SE_LSMOFA', 'MSE_MEFISTO', 'SE_MEFISTO', 'p'), 
     file = 'result/sim_S3_MSE_summary.RData')

