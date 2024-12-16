library(tidyverse)
library(refund)
library(survival)
library(tdROC)
library(ipred)
library(GPArotation)

unloadNamespace('MOFA2')


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

source('R_code/sim_S2/data_gen_MEFISTO.R')
source('source_code/functions.R')

result_MEFISTO = vector('list', 100)
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
  
  fname = paste0('result/sim_S2/fit/factor_', seed, '.RData')
  load(fname)
  
  ##
  n.factors = 3
  df_all = factors_df
  
  df_all$ID = NA
  for (i in 1:nrow(df_all)){
    tmp.sample = df_all$sample[i]
    detect = str_locate_all(tmp.sample, "_")[[1]]
    tmp.ID = str_sub(tmp.sample, detect[1, 1] + 1, detect[2, 1] - 1)
    df_all$ID[i] = as.numeric(tmp.ID)
  }
  df_all = df_all %>% arrange(ID, time)
  
  ydata_all = data.frame(.id = NA, .index = NA, .value = NA)
  ydata_all = ydata_all[-1, ]
  scores = NULL
  L1 = rep(NA, n.factors)
  for (f in 1:n.factors){
    ydata <- data.frame(.id = df_all$ID, .index = df_all$time, .value = df_all[, f])
    fpca.obj = fpca.sc(ydata = ydata, pve = 0.95, npc = 2, cov.est.method = 1)
    L1[f] = fpca.obj$npc
    scores = cbind(scores, fpca.obj$scores)
  }
  
  ### Survival estimation
  surv.dat2 = cbind(surv.dat[, c(1, 2)], scores)
  colnames(surv.dat2) = c('surv_time', 'status', 
                          paste0('scores_', 1:ncol(scores)))
  cox.obj2 = coxph(Surv(surv_time, status) ~ ., data = surv.dat2)
  p.FA = predict(cox.obj2, type = 'survival')
  
  ##
  surv.obj <- Surv(surv.dat$surv_time, surv.dat$status)
  
  starting.time = c(1:20)/20
  AUC.MEFISTO = BS.MEFISTO = rep(NA, length(starting.time))
  for (tau in starting.time){
    ## Estimation
    ROC.est <- tdROC(X = p.FA, Y = surv.dat$surv_time,
                     delta = surv.dat$status, tau = tau,
                     span = 0.1, alpha = 0.05,
                     n.grid = 1000, cut.off = 0.5)
    AUC.MEFISTO[which(tau==starting.time)] <- ROC.est$main_res$AUC.integral
    #AUC.MEFISTO[which(tau==starting.time)] <- ROC.est$AUC$value
    
    BS.MEFISTO[which(tau==starting.time)] <- sbrier(surv.obj, 1-p.FA, btime = tau)
  }
  
  result_MEFISTO[[seed]] = list(AUC.MEFISTO = AUC.MEFISTO, BS.MEFISTO = BS.MEFISTO)
  

}  ## End loop of seed 


fname = paste0('result/sim_S2/MEFISTO_predict', '.RData')
save(list = 'result_MEFISTO', file = fname)
