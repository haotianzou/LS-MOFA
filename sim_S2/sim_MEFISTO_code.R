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

source('R_code/sim_S2/data_gen_MEFISTO.R')
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
  
  metadata_MEF = data.frame(ID = NA, group = NA, sample = NA, 
                              time = NA, view = NA)[0, ]
  MEF = data.frame(group = NA, time = NA, feature = NA, 
                   value = NA, view = NA, sample = NA)[0, ]
  view = paste0('view_', 1:V)
  for (i in 1:N){
    
    Yi = time_i = vector('list', V)
    n_features = Ki = rep(NA, V)
    meta_i_all = data.frame(ID = NA, time = NA, view = NA)[0, ]
    features = vector('list', V)
    features_all = NULL
    which_time = NULL
    for (v in 1:V){
      metadata = metadata_all[[v]]
      Y = Y_all[[v]]
      
      meta_i = metadata[metadata$ID==i, ]
      meta_i$view = view[v]
      meta_i_all = rbind(meta_i_all, meta_i)
      
      time_i[[v]] = meta_i$time
      Ki[v] = nrow(meta_i)
      Yi[[v]] = matrix(Y[metadata$ID==i, ], nrow = Ki[v])
      
      n_features[v] = ncol(Y)
      features[[v]] = paste0('feature_', 1:n_features[v])
      features_all = c(features_all, rep(features[[v]], Ki[v]))
      which_time = c(which_time, 1:Ki[v])
    }
    
    ## Construct the sample pool
    time_all = unlist(time_i)
    time_all = sort(time_all[!duplicated(time_all)])
    sample_df = data.frame(time = time_all, sample = 1:length(time_all))
    
    ## Create data 
    view_i = rep(1:V, Ki)
    start = 1
    end = 0
    
    data2_metadata = data.frame(ID = meta_i_all$ID,
                                group = paste0("Subj_", meta_i_all$ID), 
                                sample = NA, 
                                time = meta_i_all$time, 
                                view = paste0('view_', view_i))
    data2 = data.frame(group = NA, time = NA, feature = features_all, 
                       value = NA, view = NA, sample = NA)
    
    for (ii in 1:sum(Ki)){
      tmp.time = data2_metadata$time[ii]
      tmp.sample = sample_df[which(sample_df$time==tmp.time), ]$sample
      data2_metadata$sample[ii] = paste0("Subj_", meta_i_all$ID[ii], "_", tmp.sample)
      
      ## 
      end = end + n_features[view_i[ii]]
      data2[start:end, ]$group = paste0("Subj_", data2_metadata$ID[ii])
      data2[start:end, ]$time = data2_metadata$time[ii]
      data2[start:end, ]$value = Yi[[view_i[ii]]][which_time[ii], ]
      data2[start:end, ]$view = paste0("view_", view_i[ii])
      data2[start:end, ]$sample = data2_metadata$sample[ii]
      
      start = end + 1
    }
    
    metadata_MEF = rbind(metadata_MEF, data2_metadata)
    MEF = rbind(MEF, data2)
  }
  
  metadata_MEF = metadata_MEF %>% arrange(view, group, time)
  MEF = MEF %>% arrange(view, group, time, feature)
  
  MOFAobject_untrained <- create_mofa(data = MEF)
  metadata_MOFA = samples_metadata(MOFAobject_untrained)
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
  
  outfile <- paste0("result/sim_S2/fit/", seed, '.hdf5')
  MOFAobject <- run_mofa(MOFAobject_untrained, outfile = outfile)
  
  factors = get_factors(MOFAobject)
  factors_df = data.frame(factor1 = NA, factor2 = NA, factor3 = NA)
  factors_df = factors_df[0, ]
  
  for (i in 1:length(factors)){
    tmp = factors[[i]]
    factors_df = rbind(factors_df, tmp)
  }
  
  factors_df$sample = rownames(factors_df)
  factors_df$time = metadata_MOFA$time
  
  fname = paste0('result/sim_S2/fit/factor_', seed, '.RData')
  save(factors_df, file = fname)
  
}  ## End loop of seed 

