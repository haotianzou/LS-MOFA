library(tidyverse)
library(MOFA2)
library(refund)
library(survival)
library(tdROC)
library(ipred)
library(GPArotation)


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

source('R_code/sim_S1/data_gen_MEFISTO.R')
source('source_code/functions.R')

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
  
  #### MEFISTO  ####
  view = paste0('view_', 1:V)
  
  meta2 = vector('list', V)
  for (v in 1:V){
    metadata = metadata_all[[v]]
    ID = unique(metadata$ID)
    
    Y = Y_all[[v]]
    nobs = nrow(metadata)
    n_features = ncol(Y)
    features = paste0('feature_', 1:n_features, '_view_', v)
    
    sample = NULL
    for (i in 1:N){
      tmp.data = metadata[metadata$ID==ID[i], ]
      sample = c(sample, c(1:nrow(tmp.data)) - 1)
    }
    
    data2_metadata = data.frame(ID = metadata$ID,
                                group = paste0("Subj_", metadata$ID), 
                                sample = paste0("Subj_", metadata$ID, "_", sample), 
                                time = metadata$time)
    meta2[[v]] = data2_metadata
    data2 = data.frame(group = NA, time = NA, feature = rep(features, nobs), 
                       value = NA, view = view[v], sample = NA)
    
    for (i in 1:nobs){
      start = (i-1)*n_features + 1
      end = i*n_features
      
      data2$group[start:end] = paste0("Subj_", metadata$ID[i])
      data2$time[start:end] = metadata$time[i]
      data2$value[start:end] = Y[i, 1:ncol(Y)]
      data2$sample[start:end] = paste0(data2$group[start:end], "_", sample[i])
      
      #if (i %% 100==0) cat(sprintf("i = %d\n", i))
    }
    
    data2$value = as.numeric(data2$value)
    
    fname = paste0('result/sim_S1/data/MEFISTO_view_', v, '_', seed, '.RData')
    save(list = 'data2', file = fname)
    
  }
  
  meta_join = full_join(meta2[[1]], meta2[[2]], by = 'sample')
  meta_join$time = ifelse(is.na(meta_join$time.x), meta_join$time.y, meta_join$time.x)
  meta_join = meta_join %>% select(sample, time)
  
  meta_join$ID = NA
  for (i in 1:nrow(meta_join)){
    tmp.sample = meta_join$sample[i]
    detect = str_locate_all(tmp.sample, "_")[[1]]
    tmp.ID = str_sub(tmp.sample, detect[1, 1] + 1, detect[2, 1] - 1)
    meta_join$ID[i] = as.numeric(tmp.ID)
  }
  
  for (v in 1:V){
    fname = paste0('result/sim_S1/data/MEFISTO_view_', v, '_', seed, '.RData')
    load(fname)
    
    if (v==1) df = data2 else df = rbind(df, data2)
  }
  
  MOFAobject_untrained <- create_mofa(data = df)
  metadata_MOFA = samples_metadata(MOFAobject_untrained) %>% left_join(meta_join, by = c('sample', 'time'))
  df2 = matrix(NA, nrow = 1, ncol = nrow(metadata_MOFA))
  df2[1, ] = metadata_MOFA$time
  colnames(df2) = metadata_MOFA$sample
  rownames(df2) = "time"
  
  MOFAobject_untrained <- set_covariates(MOFAobject_untrained,
                                         covariates = df2)
  
  ##
  data_opts <- get_default_data_options(MOFAobject_untrained)
  data_opts$center_groups <- TRUE
  
  model_opts <- get_default_model_options(MOFAobject_untrained)
  model_opts$num_factors <- 3
  n.factors = model_opts$num_factors
  
  mefisto_opts <- get_default_mefisto_options(MOFAobject_untrained)
  mefisto_opts$n_grid <- 10
  mefisto_opts$start_opt <- 20
  mefisto_opts$opt_freq <- 20
  mefisto_opts$model_groups <- FALSE # fast option (no optimization of group relationships) 
  mefisto_opts$new_values <- matrix(0:10, nrow =1) # set time points to interpolate factors to
  
  train_opts <- get_default_training_options(MOFAobject_untrained)
  train_opts$seed <- 2020
  train_opts$maxiter <- 100
  
  MOFAobject_untrained <- prepare_mofa(
    object = MOFAobject_untrained,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts,
    mefisto_options = mefisto_opts
  ) 
  
  outfile <- paste0("result/sim_S1/fit/", seed, '.hdf5')
  MOFAobject <- run_mofa(MOFAobject_untrained, outfile = outfile)
  
  # 
  factors = get_factors(MOFAobject)
  factors_df = data.frame(factor1 = NA, factor2 = NA, factor3 = NA)
  factors_df = factors_df[0, ]
  
  for (i in 1:length(factors)){
    tmp = factors[[i]]
    factors_df = rbind(factors_df, tmp)
  }
  
  factors_df$sample = rownames(factors_df)
  factors_df = factors_df %>% inner_join(metadata_MOFA, by = 'sample')
  
  fname = paste0('result/sim_S1/fit/factor_', seed, '.RData')
  save(factors_df, file = fname)
  
  ##
  
}  ## End loop of seed 

