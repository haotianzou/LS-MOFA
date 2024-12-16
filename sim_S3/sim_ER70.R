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

N = 400
V = 2
P_all = rep(200, V)
Fl = rep(2, V)
F0 = 1
L = list(L1 = rep(2, Fl[[1]] + F0), L2 = rep(2, Fl[[2]] + F0))
Lu = list(L0 = rep(2, F0))

## Define true functions
phi1 = list(phi11 = function(t) sqrt(2)*sin(pi*t), 
            phi12 = function(t) sqrt(2)*cos(pi*t))
phi2 = list(phi11 = function(t) sqrt(2)*sin(2*pi*t), 
            phi12 = function(t) sqrt(2)*cos(2*pi*t))
psi1 = list(psi1 = function(t) sqrt(2)*cos(pi*t), 
            psi2 = function(t) sqrt(2)*sin(pi*t))
phi_v1 = list(phi1 = phi1, phi2 = phi2, psi1 = psi1)

phi1 = list(phi11 = function(t) sqrt(2)*sin(pi*t), 
            phi12 = function(t) sqrt(2)*cos(pi*t))
phi2 = list(phi11 = function(t) sqrt(2)*sin(2*pi*t), 
            phi12 = function(t) sqrt(2)*cos(2*pi*t))
phi_v2 = list(phi1 = phi1, phi2 = phi2, psi1 = psi1)

phi = list(phi_v1 = phi_v1, phi_v2 = phi_v2)
psi = list(psi1 = psi1)

phi.sign_v1 <- list(phi1 = c(51, 1), phi2 = c(26, 1), psi1 = c(1, 51))
phi.sign_v2 <- list(phi1 = c(51, 1), phi2 = c(26, 1), psi1 = c(1, 51))
phi.sign = list(phi.sign_v1 = phi.sign_v1, phi.sign.v2 = phi.sign_v2)
psi.sign = list(psi1 = c(1, 51))

## Define the refined grid
t.argvals = seq(0, 1, by = 0.01)
lt = length(t.argvals)
uID_all = 1:N

source('R_code/sim_S3/data_gen_2omics_ER70.R')
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
  
  newdata = data.frame(t0 = rep(t.argvals, each = lt), t1 = rep(t.argvals, lt))
  newdata_diag = data.frame(t0 = t.argvals)
  
  #### Compute the view-wise correlation matrix C_{1p'}^(1, v')(s, t), v' = 2, ... V ####
  Y1 = Y_all[[1]]
  P1 = ncol(Y1)
  ID1 = ID_all[[1]]
  uID_1 = unique(ID1)
  time1 = time_all[[1]]
  metadata1 = metadata_all[[1]]
  Y1_base = Y1[metadata1$time==0, ]
  cov_view = vector('list', V)
  cov_base = vector('list', V)
  
  for (v in 2:V){
    Y = Y_all[[v]]
    ID = ID_all[[v]]
    uID_v = unique(ID)
    time = time_all[[v]]
    metadata = metadata_all[[v]]
    P = ncol(Y)
    nobs = table(ID)
    Y_base = Y[metadata$time==0, ]
    
    ## Compute baseline covariance ##
    cov_base[[v]] = cov(Y1_base, Y_base)
    
    cov_view[[v]] = array(NA, dim = c(P, lt, lt))
    ## Compute residual product
    Yprod_view = matrix(NA, nrow = 0, ncol = 2+P)
    colnames(Yprod_view) = c('t0', 't1', paste0('Yprod_1_', 1:P))
    
    for (i in 1:N){
      tmp.Y1 = matrix(Y1[which(metadata1$ID==uID_1[i]), ], ncol = P1)
      tmp.time1 = time1[which(metadata1$ID==uID_1[i])]
      tmp.lt1 = length(tmp.time1)
      
      tmp.Y = matrix(Y[which(metadata$ID==uID_v[i]), ], ncol = P)
      tmp.time = time[which(metadata$ID==uID_v[i])]
      tmp.lt = length(tmp.time)
      
      t0 = rep(tmp.time1, each = tmp.lt)
      t1 = rep(tmp.time, tmp.lt1)
      
      tmp.data = matrix(NA, nrow = tmp.lt1*tmp.lt, ncol = 2+P)
      colnames(tmp.data) = c('t0', 't1', paste0('Yprod_', 1:P))
      tmp.data[, 1] = t0
      tmp.data[, 2] = t1
      
      for (p in 1:P){
        tmp.data[, p+2] = rep(tmp.Y1[, 1], each = tmp.lt)*rep(tmp.Y[, p], tmp.lt1)
      }
      
      Yprod_view = rbind(Yprod_view, tmp.data)
    }
    
    Yprod_view = as.data.frame(Yprod_view)
    
    ## view-wise covariance estimation
    for (p in 1:P){
      tmp.Yprod = Yprod_view[, c(1, 2, 2 + p)]
      colnames(tmp.Yprod)[3] = 'Yprod'
      
      gam.obj = gam(Yprod ~ te(t0, t1, bs = 'ts'), data = tmp.Yprod)
      pred = predict(gam.obj, newdata)
      tmp.cov = matrix(pred, nrow = lt, ncol = lt, byrow = T)
      tmp.cov = (tmp.cov + t(tmp.cov))/2
      cov_view[[v]][p, , ] = tmp.cov
    }
  }
  
  Lambda.est = vector('list', V)
  for (v in 2:V){
    tmp.cov = cov_view[[v]][, 1, 1]
    Lambda.est[[v]] = matrix(tmp.cov, ncol = 1)
  }
  
  #### Compute W and C for v = 2, ..., V ####
  sigma_p_all = vector('list', V)
  cov_est_result = vector('list', V)
  for (v in 2:V) {
    time = time_all[[v]]
    ID = ID_all[[v]]
    
    uID = unique(ID)
    nobs = table(ID)
    k0 = unname(c(0, cumsum(nobs)[1:(N-1)]) + 1)
    k1 = unname(cumsum(nobs))
    
    ## 
    Y = Y_all[[v]]
    metadata = metadata_all[[v]]
    P = ncol(Y)
    
    ## Baseline information
    Y_base = Y[metadata$time==0, ]
    cov.mat = cov(Y_base)
    cov.mat = cov.mat - Lambda.est[[v]] %*% t(Lambda.est[[v]])
    
    cov.mat.diag = matrix(0, nrow = P, ncol = P)
    diag(cov.mat.diag) = diag(cov.mat)
    
    ## Start with an initial matrix of Psi ##
    Psi = matrix(0, nrow = P, ncol = P)
    diag(Psi) = rep(0.01, P)
    
    ## Iteration ##
    max.iters = 200
    curr.iter = 0
    convergence = F
    eps = 1e-6
    while (!convergence){
      curr.iter = curr.iter + 1
      W.prod.curr = cov.mat - Psi
      eigen.obj = eigen(W.prod.curr)
      values = eigen.obj$values[1:Fl[v]]
      vectors = eigen.obj$vectors[, 1:Fl[v]]
      
      values_sqrt = matrix(0, nrow = Fl[v], ncol = Fl[v])
      diag(values_sqrt) = sqrt(values) 
      
      W.curr = vectors %*% values_sqrt
      W.curr = Varimax(W.curr)$loadings
      W.curr = abs(W.curr)
      
      W.prod.curr = W.curr %*% t(W.curr)
      W.prod.curr.diag = matrix(0, nrow = P, ncol = P)
      diag(W.prod.curr.diag) = diag(W.prod.curr)
      
      Psi.curr = cov.mat.diag - W.prod.curr.diag
      ## Psi to be capped at 0.01 ##
      diag.Psi.curr = diag(Psi.curr)
      diag.Psi.curr = ifelse(diag.Psi.curr<0, 0.01, diag.Psi.curr)
      diag(Psi.curr) = diag.Psi.curr
      
      diff = mean((diag(Psi.curr) - diag(Psi))^2)
      
      if (abs(diff)<eps) {
        convergence = T
        cat(sprintf('After %d iterations, Model converged!\n', curr.iter))
        break;
      } 
      
      Psi = Psi.curr
    }
    
    sigma_p_all[[v]] = sqrt(diag(Psi))
    
    ## Covariance estimation
    Yprod = matrix(NA, nrow = 0, ncol = 2+P)
    Yprod_cross = matrix(NA, nrow = 0, ncol = 2+P)
    colnames(Yprod) = c('t0', 't1', paste0('Yprod_', 1:P))
    colnames(Yprod_cross) = c('t0', 't1', paste0('Yprod_1_', 1:P))
    
    for (i in 1:N){
      tmp.Y = matrix(Y[which(metadata$ID==uID[i]), ], ncol = P)
      tmp.time = time[which(metadata$ID==uID[i])]
      tmp.lt = length(tmp.time)
      t0 = rep(tmp.time, each = tmp.lt)
      t1 = rep(tmp.time, tmp.lt)
      
      tmp.data = matrix(NA, nrow = tmp.lt*tmp.lt, ncol = 2+P)
      colnames(tmp.data) = c('t0', 't1', paste0('Yprod_', 1:P))
      tmp.data[, 1] = t0
      tmp.data[, 2] = t1
      
      for (p in 1:P){
        tmp.data[, p+2] = rep(tmp.Y[, p], each = tmp.lt)*rep(tmp.Y[, p], tmp.lt)
      }
      
      Yprod = rbind(Yprod, tmp.data)
      
      tmp.data = matrix(NA, nrow = tmp.lt*tmp.lt, ncol = 2+P)
      colnames(tmp.data) = c('t0', 't1', paste0('Yprod_1_', 1:P))
      tmp.data[, 1] = t0
      tmp.data[, 2] = t1
      
      for (p in 1:P){
        tmp.data[, p+2] = rep(tmp.Y[, 1], each = tmp.lt)*rep(tmp.Y[, p], tmp.lt)
      }
      
      Yprod_cross = rbind(Yprod_cross, tmp.data)
    }
    
    Yprod = as.data.frame(Yprod)
    Yprod_cross = as.data.frame(Yprod_cross)
    
    Yprod_diag = Yprod[which(Yprod[, 1]==Yprod[, 2]), ]
    Yprod_offdiag = Yprod[which(Yprod[, 1]!=Yprod[, 2]), ] ## extract off diagonal elements
    
    ## Auto- and cross-covariance
    cov_auto = array(NA, dim = c(P, lt, lt))
    cov_cross = array(NA, dim = c(P, lt, lt))
    for (p in 1:P){
      ## Auto-covariance
      tmp.Yprod_offdiag = Yprod_offdiag[, c(1, 2, 2 + p)]
      colnames(tmp.Yprod_offdiag)[3] = 'Yprod'
      
      gam.obj = gam(Yprod ~ te(t0, t1, bs = 'ts'), data = tmp.Yprod_offdiag)
      pred = predict(gam.obj, newdata)
      tmp.cov = matrix(pred, nrow = lt, ncol = lt, byrow = T)
      tmp.cov = (tmp.cov + t(tmp.cov))/2
      cov_auto[p, , ] = tmp.cov
      
      ## Cross-covariance
      tmp.Yprod_cross = Yprod_cross[, c(1, 2, 2 + p)]
      colnames(tmp.Yprod_cross)[3] = 'Yprod'
      
      gam.obj = gam(Yprod ~ te(t0, t1, bs = 'ts'), data = tmp.Yprod_cross)
      pred = predict(gam.obj, newdata)
      tmp.cov = matrix(pred, nrow = lt, ncol = lt, byrow = T)
      tmp.cov = (tmp.cov + t(tmp.cov))/2
      cov_cross[p, , ] = tmp.cov
    }
    
    ## Estimation of covariance matrix 
    C_mat = array(NA, dim = c(Fl[v], lt, lt))
    for (t0 in 1:lt){
      for (t1 in 1:lt){
        tmp.cov1 = cov_auto[, t0, t1] - Lambda.est[[v]]*cov_view[[v]][, t0, t1]
        tmp.cov2 = cov_cross[2:P, t0, t1] - Lambda.est[[v]][1]*cov_view[[v]][2:P, t0, t1]
        
        tmp.Cy = matrix(c(tmp.cov1, tmp.cov2), ncol = 1)
        x3 = cbind(1, rbind(W.curr[, ]^2, W.curr[1, ]*W.curr[2:P, ]))
        x3 = ifelse(x3<0.01, 0, x3)
        coef.lm = solve(t(x3) %*% x3, t(x3)) %*% tmp.Cy
        C_mat[1:Fl[v], t0, t1] = coef.lm[2:(1+Fl[v])]
      }
    }
    
    
    ## 
    first.eigenvalue = rep(NA, Fl[v])
    for (f in 1:Fl[v]){
      eigen.obj = eigen(C_mat[f, , ])
      first.eigenvalue[f] = eigen.obj$values[1]/lt
    }
    
    ## Swap the covariance matrix so that the first eigenvalue is ordered from largest to smallest
    for (f0 in 1:(Fl[v] - 1)){
      for (f1 in (f0 + 1):Fl[v]){
        if (first.eigenvalue[f0]<first.eigenvalue[f1]){
          tmp = C_mat[f0, , ]
          C_mat[f0, , ] = C_mat[f1, , ]
          C_mat[f1, , ] = tmp
          
          tmp = W.curr[, f0]
          W.curr[, f0] = W.curr[, f1]
          W.curr[, f1] = tmp
        }
      }
    }
    
    eigendiff_within = rep(NA, Fl[v])
    eigendiff_cross = rep(NA, Fl[v] - 1)
    first.eigenvalue = rep(NA, Fl[v])
    for (f in 1:Fl[v]){
      eigen.obj = eigen(C_mat[f, , ])
      values = eigen.obj$values[1:L[[v]][f]]/lt
      first.eigenvalue[f] = values[1]
      eigendiff_within[f] = values[1] - values[2]
    }
    
    for (f in 1:(Fl[v] - 1)){
      eigendiff_cross[f] = first.eigenvalue[f] - first.eigenvalue[f+1]
    }
    
    eigendiff = c(eigendiff_cross, eigendiff_within)
    if (sum(eigendiff<0.5)>0) convergence = F; 
    
    phi.sign_v = phi.sign[[v]]
    phi_v = phi[[v]]
    
    ## Estimate eigenvalues and eigenfunctions
    phi_est = vector('list', Fl[v])
    d_est = vector('list', Fl[v])
    MSE = vector('list', Fl[v])
    for (f in 1:Fl[v]){
      eigen.obj = eigen(C_mat[f, , ])
      
      values = eigen.obj$values[eigen.obj$values>0]
      # pve = cumsum(values)/sum(values)
      # L[f] = min(which(pve>0.95))
      L0 = L[[v]][f]
      
      values = values[1:L0]/lt
      vectors = matrix(eigen.obj$vectors[, 1:L0]*sqrt(lt), ncol = L0)
      
      phi_est[[f]] = matrix(NA, nrow = lt, ncol = L0)
      d_est[[f]] = values
      MSE[[f]] = rep(NA, L0)
      for (l in 1:L0){
        if (sign(vectors[phi.sign_v[[f]][l], l])!=1) vectors[, l] = -vectors[, l]
        phi_est[[f]][, l] = vectors[, l]
        MSE[[f]][l] = mean((phi_est[[f]][, l] - phi_v[[f]][[l]](t.argvals))^2)
      }
    }
    
    cov_est_result[[v]] = list(phi_est = phi_est, d_est = d_est, MSE = MSE, W_est = W.curr)
  }
  
  if (!convergence) next; 
  
  #### Estimate C0(s, t) ####
  t0 = rep(t.argvals, each = lt)
  t1 = rep(t.argvals, lt)
  df_all = data.frame(t0 = NA, t1 = NA, value = NA)
  df_all = df_all[0, ]
  
  for (v in 2:V){
    Y = Y_all[[v]]
    metadata = metadata_all[[v]]
    P = ncol(Y)
    for (p in 1:P){
      if (abs(Lambda.est[[v]][p])>0.5){
        tmp = cov_view[[v]][p, , ]/Lambda.est[[v]][p]
        df = data.frame(t0 = t0, t1 = t1, value = as.vector(tmp))
        df_all = rbind(df_all, df)
      }
    }
  }
  
  if (nrow(df_all)/(lt*lt)<=5){
    df_all = data.frame(t0 = NA, t1 = NA, value = NA)
    df_all = df_all[0, ]
    
    for (v in 2:V){
      Y = Y_all[[v]]
      metadata = metadata_all[[v]]
      P = ncol(Y)
      for (p in 1:P){
        if (abs(Lambda.est[[v]][p])>0.3){
          tmp = cov_view[[v]][p, , ]/Lambda.est[[v]][p]
          df = data.frame(t0 = t0, t1 = t1, value = as.vector(tmp))
          df_all = rbind(df_all, df)
        }
      }
    }
    
  }
  
  gam.obj = gam(value ~ te(t0, t1, bs = 'ts'), data = df_all)
  pred = predict(gam.obj, newdata)
  tmp.cov = matrix(pred, nrow = lt, ncol = lt, byrow = T)
  tmp.cov = (tmp.cov + t(tmp.cov))/2
  C0_est_bak = tmp.cov
  
  ## To verify
  # C0_est = C[[1]][[3]]
  C0_est = C0_est_bak
  # mean(diff^2)
  
  #### Compute MSE of psi ####
  psi_est = vector('list', F0)
  d0_est = vector('list', F0)
  MSE.psi = vector('list', F0)
  for (f in 1:F0){
    eigen.obj = eigen(C0_est)
    
    values = eigen.obj$values[eigen.obj$values>0]
    # pve = cumsum(values)/sum(values)
    # L[f] = min(which(pve>0.95))
    
    values = values[1:Lu[[f]]]/lt
    vectors = matrix(eigen.obj$vectors[, 1:Lu[[f]]]*sqrt(lt), ncol = Lu[[f]])
    
    psi_est[[f]] = matrix(NA, nrow = lt, ncol = Lu[[f]])
    d0_est[[f]] = values
    MSE.psi[[f]] = rep(NA, Lu[[f]])
    for (l in 1:Lu[[f]]){
      if (sign(vectors[psi.sign[[f]][l], l])!=1) vectors[, l] = -vectors[, l]
      psi_est[[f]][, l] = vectors[, l]
      if (!is.null(psi)){
        MSE.psi[[f]][l] = mean((psi_est[[f]][, l] - psi[[f]][[l]](t.argvals))^2)
      }
    }
    
    d0_mat = matrix(0, nrow = Lu[[f]], ncol = Lu[[f]])
    diag(d0_mat) = d0_est[[f]]
    
    C0_est = psi_est[[f]] %*% d0_mat %*% t(psi_est[[f]])
  }
  
  #### For v = 1, estimate covariance matrix ####
  for (v in 1:1) {
    time = time_all[[v]]
    ID = ID_all[[v]]
    
    uID = unique(ID)
    nobs = table(ID)
    k0 = unname(c(0, cumsum(nobs)[1:(N-1)]) + 1)
    k1 = unname(cumsum(nobs))
    
    ## 
    Y = Y_all[[v]]
    metadata = metadata_all[[v]]
    P = ncol(Y)
    
    ## for v = 2, ..., V, estimate lambda_{1p} ##
    Lambda.est[[v]] = matrix(NA, nrow = P, ncol = F0)
    Lambda.est[[v]][1, ] = 1
    for (p in 2:P){
      tmp.Y = NULL
      tmp.X = NULL
      for (v0 in 2:V){
        tmp.Y = c(tmp.Y, cov_base[[v0]][p, ])
        tmp.X = c(tmp.X, Lambda.est[[v0]])
      }
      tmp.X = cbind(1, tmp.X)
      coef.lm = solve(t(tmp.X) %*% tmp.X, t(tmp.X)) %*% tmp.Y
      Lambda.est[[v]][p, ] = coef.lm[2]
    }
    
    diff = Lambda.est[[v]] - W[[v]][, 3]
    
    ## Estimation of W 
    ## Baseline information
    Y_base = Y[metadata$time==0, ]
    cov.mat = cov(Y_base)
    cov.mat = cov.mat - Lambda.est[[v]] %*% t(Lambda.est[[v]])
    
    cov.mat.diag = matrix(0, nrow = P, ncol = P)
    diag(cov.mat.diag) = diag(cov.mat)
    
    ## Start with an initial matrix of Psi ##
    Psi = matrix(0, nrow = P, ncol = P)
    diag(Psi) = rep(0.01, P)
    
    ## Iteration ##
    max.iters = 200
    curr.iter = 0
    convergence = F
    eps = 1e-6
    while (!convergence){
      curr.iter = curr.iter + 1
      W.prod.curr = cov.mat - Psi
      eigen.obj = eigen(W.prod.curr)
      values = eigen.obj$values[1:Fl[v]]
      vectors = eigen.obj$vectors[, 1:Fl[v]]
      
      values_sqrt = matrix(0, nrow = Fl[v], ncol = Fl[v])
      diag(values_sqrt) = sqrt(values) 
      
      W.curr = vectors %*% values_sqrt
      W.curr = Varimax(W.curr)$loadings
      W.curr = abs(W.curr)
      
      W.prod.curr = W.curr %*% t(W.curr)
      W.prod.curr.diag = matrix(0, nrow = P, ncol = P)
      diag(W.prod.curr.diag) = diag(W.prod.curr)
      
      Psi.curr = cov.mat.diag - W.prod.curr.diag
      ## Psi to be capped at 0.01 ##
      diag.Psi.curr = diag(Psi.curr)
      diag.Psi.curr = ifelse(diag.Psi.curr<0, 0.01, diag.Psi.curr)
      diag(Psi.curr) = diag.Psi.curr
      
      diff = mean((diag(Psi.curr) - diag(Psi))^2)
      
      if (abs(diff)<eps) {
        convergence = T
        cat(sprintf('After %d iterations, Model converged!\n', curr.iter))
        break;
      } 
      
      Psi = Psi.curr
    }
    
    if (!convergence) next; 
    sigma_p_all[[v]] = sqrt(diag(Psi))
    
    ## Covariance estimation
    Yprod = matrix(NA, nrow = 0, ncol = 2+P)
    Yprod_cross = matrix(NA, nrow = 0, ncol = 2+P)
    colnames(Yprod) = c('t0', 't1', paste0('Yprod_', 1:P))
    colnames(Yprod_cross) = c('t0', 't1', paste0('Yprod_1_', 1:P))
    
    for (i in 1:N){
      tmp.Y = matrix(Y[which(metadata$ID==uID[i]), ], ncol = P)
      tmp.time = time[which(metadata$ID==uID[i])]
      tmp.lt = length(tmp.time)
      t0 = rep(tmp.time, each = tmp.lt)
      t1 = rep(tmp.time, tmp.lt)
      
      tmp.data = matrix(NA, nrow = tmp.lt*tmp.lt, ncol = 2+P)
      colnames(tmp.data) = c('t0', 't1', paste0('Yprod_', 1:P))
      tmp.data[, 1] = t0
      tmp.data[, 2] = t1
      
      for (p in 1:P){
        tmp.data[, p+2] = rep(tmp.Y[, p], each = tmp.lt)*rep(tmp.Y[, p], tmp.lt)
      }
      
      Yprod = rbind(Yprod, tmp.data)
      
      tmp.data = matrix(NA, nrow = tmp.lt*tmp.lt, ncol = 2+P)
      colnames(tmp.data) = c('t0', 't1', paste0('Yprod_1_', 1:P))
      tmp.data[, 1] = t0
      tmp.data[, 2] = t1
      
      for (p in 1:P){
        tmp.data[, p+2] = rep(tmp.Y[, 1], each = tmp.lt)*rep(tmp.Y[, p], tmp.lt)
      }
      
      Yprod_cross = rbind(Yprod_cross, tmp.data)
    }
    
    Yprod = as.data.frame(Yprod)
    Yprod_cross = as.data.frame(Yprod_cross)
    
    Yprod_diag = Yprod[which(Yprod[, 1]==Yprod[, 2]), ]
    Yprod_offdiag = Yprod[which(Yprod[, 1]!=Yprod[, 2]), ] ## extract off diagonal elements
    
    ## Auto- and cross-covariance
    cov_auto = array(NA, dim = c(P, lt, lt))
    cov_cross = array(NA, dim = c(P, lt, lt))
    for (p in 1:P){
      ## Auto-covariance
      tmp.Yprod_offdiag = Yprod_offdiag[, c(1, 2, 2 + p)]
      colnames(tmp.Yprod_offdiag)[3] = 'Yprod'
      
      gam.obj = gam(Yprod ~ te(t0, t1, bs = 'ts'), data = tmp.Yprod_offdiag)
      pred = predict(gam.obj, newdata)
      tmp.cov = matrix(pred, nrow = lt, ncol = lt, byrow = T)
      tmp.cov = (tmp.cov + t(tmp.cov))/2
      cov_auto[p, , ] = tmp.cov
      
      ## Cross-covariance
      tmp.Yprod_cross = Yprod_cross[, c(1, 2, 2 + p)]
      colnames(tmp.Yprod_cross)[3] = 'Yprod'
      
      gam.obj = gam(Yprod ~ te(t0, t1, bs = 'ts'), data = tmp.Yprod_cross)
      pred = predict(gam.obj, newdata)
      tmp.cov = matrix(pred, nrow = lt, ncol = lt, byrow = T)
      tmp.cov = (tmp.cov + t(tmp.cov))/2
      cov_cross[p, , ] = tmp.cov
    }
    
    
    ## Estimation of covariance matrix 
    C_mat = array(NA, dim = c(Fl[v], lt, lt))
    for (t0 in 1:lt){
      for (t1 in 1:lt){
        tmp.cov1 = cov_auto[, t0, t1] - Lambda.est[[v]]^2*C0_est[t0, t1]
        tmp.cov2 = cov_cross[2:P, t0, t1] - Lambda.est[[v]][1]*Lambda.est[[v]][2:P]*C0_est[t0, t1]
        
        tmp.Cy = matrix(c(tmp.cov1, tmp.cov2), ncol = 1)
        x3 = cbind(1, rbind(W.curr[, ]^2, W.curr[1, ]*W.curr[2:P, ]))
        x3 = ifelse(x3<0.01, 0, x3)
        coef.lm = solve(t(x3) %*% x3, t(x3)) %*% tmp.Cy
        C_mat[1:Fl[v], t0, t1] = coef.lm[2:(1+Fl[v])]
      }
    }
    
    ## 
    first.eigenvalue = rep(NA, Fl[v])
    for (f in 1:Fl[v]){
      eigen.obj = eigen(C_mat[f, , ])
      first.eigenvalue[f] = eigen.obj$values[1]/lt
    }
    
    ## Swap the covariance matrix so that the first eigenvalue is ordered from largest to smallest
    for (f0 in 1:(Fl[v] - 1)){
      for (f1 in (f0 + 1):Fl[v]){
        if (first.eigenvalue[f0]<first.eigenvalue[f1]){
          tmp = C_mat[f0, , ]
          C_mat[f0, , ] = C_mat[f1, , ]
          C_mat[f1, , ] = tmp
          
          tmp = W.curr[, f0]
          W.curr[, f0] = W.curr[, f1]
          W.curr[, f1] = tmp
        }
      }
    }
    
    eigendiff_within = rep(NA, Fl[v])
    eigendiff_cross = rep(NA, Fl[v] - 1)
    first.eigenvalue = rep(NA, Fl[v])
    for (f in 1:Fl[v]){
      eigen.obj = eigen(C_mat[f, , ])
      values = eigen.obj$values[1:L[[v]][f]]/lt
      first.eigenvalue[f] = values[1]
      eigendiff_within[f] = values[1] - values[2]
    }
    
    for (f in 1:(Fl[v] - 1)){
      eigendiff_cross[f] = first.eigenvalue[f] - first.eigenvalue[f+1]
    }
    
    eigendiff = c(eigendiff_cross, eigendiff_within)
    if (sum(eigendiff<0.5)>0) convergence = F; 
    
    phi.sign_v = phi.sign[[v]]
    phi_v = phi[[v]]
    
    ## Estimate eigenvalues and eigenfunctions
    phi_est = vector('list', Fl[v])
    d_est = vector('list', Fl[v])
    MSE = vector('list', Fl[v])
    for (f in 1:Fl[v]){
      eigen.obj = eigen(C_mat[f, , ])
      
      values = eigen.obj$values[eigen.obj$values>0]
      # pve = cumsum(values)/sum(values)
      # L[f] = min(which(pve>0.95))
      L0 = L[[v]][f]
      
      values = values[1:L0]/lt
      vectors = matrix(eigen.obj$vectors[, 1:L0]*sqrt(lt), ncol = L0)
      
      phi_est[[f]] = matrix(NA, nrow = lt, ncol = L0)
      d_est[[f]] = values
      MSE[[f]] = rep(NA, L0)
      for (l in 1:L0){
        if (sign(vectors[phi.sign_v[[f]][l], l])!=1) vectors[, l] = -vectors[, l]
        phi_est[[f]][, l] = vectors[, l]
        MSE[[f]][l] = mean((phi_est[[f]][, l] - phi_v[[f]][[l]](t.argvals))^2)
      }
    }
    cov_est_result[[v]] = list(phi_est = phi_est, d_est = d_est, MSE = MSE, W_est = W.curr)
  }
  
  if (!convergence) next; 
  
  C0_result = list(psi_est = psi_est, d0_est = d0_est, MSE.psi = MSE.psi)
  
  result_all = list(cov_est_result = cov_est_result, C0_result = C0_result, Lambda.est = Lambda.est,
                    W_true = data$W)
  
  fname = paste0('result/sim_S3/LS-MOFA2/', N, '_ER70/', seed, '.RData')
  save(list = 'result_all', file = fname)
  
  cat(sprintf("seed = %d\n", seed))
  
  
}  ## End loop of seed 

