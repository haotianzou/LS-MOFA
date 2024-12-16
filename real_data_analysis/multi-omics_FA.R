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

source('source_code/functions.R')

load('dataset/new/multi-omics_imputed.RData')
V = data_all$V
N = data_all$N
Y_all = data_all$Y
metadata_all = data_all$metadata
P_all = c(ncol(Y_all[[1]]), ncol(Y_all[[2]]))  ## 749, 104

seed = 1
set.seed(2024*seed)

##
t.argvals = seq(0, 1, by = 0.01)
lt = length(t.argvals)
uID_all = 1:N

## De-mean the outcome
for (v in 1:V){
  Y = Y_all[[v]]
  time = metadata_all[[v]]$time
  P = ncol(Y)
  
  for (p in 1:P){
    df = data.frame(time = time, value = Y[, p])
    gam.obj = gam(value ~ te(time, bs = 'cs'), data = df)
    pred = predict(gam.obj)
    Y_all[[v]][, p] = Y_all[[v]][, p] - pred
  }
}

## Determine the number of factors
F0 = 1
Fl = rep(NA, V)
v = 1

## First omics 
Y_base = Y_all[[v]][metadata_all[[v]]$time==0, ]
# scree(Y_base)
Fl[v] = 2  ## 40.3% total variation
# fa.obj = fa(Y_base, nfactors = Fl[v] + F0, rotate = 'varimax')
# loadings = matrix(fa.obj$loadings, ncol = Fl[v] + F0)

## Second omics
v = 2
Y_base = Y_all[[v]][metadata_all[[v]]$time==0, ]
#scree(Y_base)
Fl[v] = 6  # 44.5% total variation
# fa.obj = fa(Y_base, nfactors = Fl[v] + F0, rotate = 'varimax')
# loadings = matrix(fa.obj$loadings, ncol = Fl[v] + F0)

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
      tmp.cov1 = cov_auto[, t0, t1] - Lambda.est[[v]]*cov_view[[v]][, t0, t1]
      tmp.cov2 = cov_cross[2:P, t0, t1] - Lambda.est[[v]][1]*cov_view[[v]][2:P, t0, t1]
      # tmp.cov2 = cov_cross[2:P, t0, t1] - Lambda.est[[v]][2:P]*cov_view[[v]][1, t0, t1]
      tmp.Cy = matrix(c(tmp.cov1, tmp.cov2), ncol = 1)
      x3 = cbind(1, rbind(W.curr[, ]^2, W.curr[1, ]*W.curr[2:P, ]))
      x3 = ifelse(x3<0.01, 0, x3)
      coef.lm = solve(t(x3) %*% x3, t(x3)) %*% tmp.Cy
      C_mat[1:Fl[v], t0, t1] = coef.lm[2:(1+Fl[v])]
    }
  }
  
  ## Order the eigenvalue from largest to smallest
  ordered.eigenvalue = FALSE
  while (!ordered.eigenvalue){
    
    ##
    first.eigenvalue = rep(NA, Fl[v])
    for (f in 1:Fl[v]){
      eigen.obj = eigen(C_mat[f, , ])
      first.eigenvalue[f] = eigen.obj$values[1]/lt
    }
    
    ordered.eigenvalue = TRUE
    for (f in 1:(Fl[v]-1)){
      if (first.eigenvalue[f]<first.eigenvalue[f+1]) ordered.eigenvalue = FALSE
    }
    
    if (ordered.eigenvalue) break; 
    
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
    
    ##
  }
  
  
  ## Estimate eigenvalues and eigenfunctions
  phi_est = vector('list', Fl[v])
  d_est = vector('list', Fl[v])
  L = rep(NA, Fl[v])
  for (f in 1:Fl[v]){
    eigen.obj = eigen(C_mat[f, , ])
    
    values = eigen.obj$values[eigen.obj$values>0]
    pve = cumsum(values)/sum(values)
    L[f] = min(which(pve>0.95))
    
    L1 = L[f]
    
    values = values[1:L1]/lt
    vectors = matrix(eigen.obj$vectors[, 1:L1]*sqrt(lt), ncol = L1)
    
    phi_est[[f]] = matrix(NA, nrow = lt, ncol = L1)
    d_est[[f]] = values
    for (l in 1:L1){
      if (sign(vectors[1, l])!=1) vectors[, l] = -vectors[, l]
      phi_est[[f]][, l] = vectors[, l]
    }
  }
  
  cov_est_result[[v]] = list(phi_est = phi_est, d_est = d_est, W_est = W.curr, L = L)
}

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
    if (abs(Lambda.est[[v]][p])>0.1){
      tmp = cov_view[[v]][p, , ]/Lambda.est[[v]][p]
      df = data.frame(t0 = t0, t1 = t1, value = as.vector(tmp))
      df_all = rbind(df_all, df)
    }
  }
}

gam.obj = gam(value ~ te(t0, t1, bs = 'ts'), data = df_all)
pred = predict(gam.obj, newdata)
tmp.cov = matrix(pred, nrow = lt, ncol = lt, byrow = T)
tmp.cov = (tmp.cov + t(tmp.cov))/2
C0_est_bak = tmp.cov
C0_est = C0_est_bak

#### Correct sign of psi ####
psi_est = vector('list', F0)
d0_est = vector('list', F0)
L0 = NULL
for (f in 1:F0){
  eigen.obj = eigen(C0_est)
  
  values = eigen.obj$values[eigen.obj$values>0]
  pve = cumsum(values)/sum(values)
  L0 = min(which(pve>0.95))
  
  values = values[1:L0]/lt
  vectors = matrix(eigen.obj$vectors[, 1:L0]*sqrt(lt), ncol = L0)
  
  psi_est[[f]] = matrix(NA, nrow = lt, ncol = L0)
  d0_est[[f]] = values
  for (l in 1:L0){
    if (sign(vectors[1, l])!=1) vectors[, l] = -vectors[, l]
    psi_est[[f]][, l] = vectors[, l]
  }
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
  
  
  ## Order the eigenvalue from largest to smallest
  ordered.eigenvalue = FALSE
  while (!ordered.eigenvalue){
    
    ##
    first.eigenvalue = rep(NA, Fl[v])
    for (f in 1:Fl[v]){
      eigen.obj = eigen(C_mat[f, , ])
      first.eigenvalue[f] = eigen.obj$values[1]/lt
    }
    
    ordered.eigenvalue = TRUE
    for (f in 1:(Fl[v]-1)){
      if (first.eigenvalue[f]<first.eigenvalue[f+1]) ordered.eigenvalue = FALSE
    }
    
    if (ordered.eigenvalue) break; 
    
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
    
    ##
  }
  
  
  ## Estimate eigenvalues and eigenfunctions
  phi_est = vector('list', Fl[v])
  d_est = vector('list', Fl[v])
  L = rep(NA, Fl[v])
  for (f in 1:Fl[v]){
    eigen.obj = eigen(C_mat[f, , ])
    
    values = eigen.obj$values[eigen.obj$values>0]
    pve = cumsum(values)/sum(values)
    L[f] = min(which(pve>0.95))
    
    L1 = L[f]
    
    values = values[1:L1]/lt
    vectors = matrix(eigen.obj$vectors[, 1:L1]*sqrt(lt), ncol = L1)
    
    phi_est[[f]] = matrix(NA, nrow = lt, ncol = L1)
    d_est[[f]] = values
    for (l in 1:L1){
      if (sign(vectors[1, l])!=1) vectors[, l] = -vectors[, l]
      phi_est[[f]][, l] = vectors[, l]
    }
  }
  
  cov_est_result[[v]] = list(phi_est = phi_est, d_est = d_est, W_est = W.curr, L = L)
}

C0_result = list(psi_est = psi_est, d0_est = d0_est, L0 = L0)

result_all = list(cov_est_result = cov_est_result, C0_result = C0_result, Lambda.est = Lambda.est)
fname = paste0('RData/new/', 'cov_est_new.RData')
save(list = 'result_all', file = fname)

#### Estimate Ui(t) and zeta ####
zeta_est = matrix(NA, nrow = N, ncol = L0)
U_est = vector('list', V)
L = vector('list', V)

for (v in 1:V){
  U_est[[v]] = vector('list', N)
  L1 = cov_est_result[[v]]$L
  L[[v]] = c(L1, L0)
}

for (i in 1:N){
  
  ## Extract essential information
  Yi_all = NULL
  P = rep(NA, V)
  Ki = rep(NA, V)
  D = vector('list', V)
  phi_obs = vector('list', V)
  time_i = vector('list', V)
  W_est = vector('list', V)
  sigma2_p = vector('list', V)
  
  ##
  for (v in 1:V){
    W_est[[v]] = cbind(cov_est_result[[v]]$W_est, Lambda.est[[v]])
    sigma2_p[[v]] = sigma_p_all[[v]]^2
    
    Y = Y_all[[v]]
    P[v] = ncol(Y)
    metadata = metadata_all[[v]]
    
    Yi = matrix(Y[metadata$ID==i, ], ncol = P[v])
    time_i[[v]] = metadata[metadata$ID==i, ]$time
    Ki[v] = nrow(Yi)
    Yi = as.vector(Yi)  ## vectorize Y: 1st column; 2nd column, etc.
    Yi_all = c(Yi_all, Yi)
    
    D[[v]] = vector('list', Fl[v] + F0)
    phi_obs[[v]] = vector('list', Fl[v] + F0)
    
    for (f in 1:Fl[v]){
      vectors = cov_est_result[[v]]$phi_est[[f]]
      values = cov_est_result[[v]]$d_est[[f]]
      D[[v]][[f]] = matrix(0, nrow = L[[v]][f], ncol = L[[v]][f])
      diag(D[[v]][[f]]) = values
      
      ## Smooth the eigenfunction on the time of each subject
      phi_obs[[v]][[f]] = matrix(NA, nrow = Ki[v], ncol = L[[v]][f])
      for (l in 1:L[[v]][f]){
        phi_obs[[v]][[f]][, l] = bs.smooth(vectors[, l], t.argvals, time_i[[v]])$est.value
      }
    }
    
    for (f in (Fl[v] + 1):(Fl[v] + F0)){
      vectors = C0_result$psi_est[[f - Fl[v]]]
      values = C0_result$d0_est[[f - Fl[v]]]
      D[[v]][[f]] = matrix(0, nrow = L[[v]][f], ncol = L[[v]][f])
      diag(D[[v]][[f]]) = values
      
      ## Smooth the eigenfunction on the time of each subject
      phi_obs[[v]][[f]] = matrix(NA, nrow = Ki[v], ncol = L[[v]][f])
      for (l in 1:L[[v]][f]){
        phi_obs[[v]][[f]][, l] = bs.smooth(vectors[, l], t.argvals, time_i[[v]])$est.value
      }
    }
    
  } # End loop of v
  
  Yi_all = matrix(Yi_all, ncol = 1)
  
  ## Covariance of Y: Sigma_11
  Sigma11 = matrix(NA, nrow = sum(Ki*P), ncol = sum(Ki*P))
  k0 = c(0, cumsum(Ki*P)[1:(V-1)]) + 1
  k1 = cumsum(Ki*P)
  
  ## For diagonal elements: v0 = v1 = v
  for (v in 1:V){
    
    tmp.Sigma11 = matrix(NA, nrow = Ki[v]*P[v], ncol = Ki[v]*P[v])
    
    tmp.C = vector('list', Fl[v] + F0)
    for (f in 1:(Fl[v] + F0)){
      phi.curr = matrix(phi_obs[[v]][[f]], nrow = Ki[v])
      tmp.C[[f]] = phi.curr %*% D[[v]][[f]] %*% t(phi.curr)
    }
    
    for (p0 in 1:P[v]){
      for (p1 in 1:P[v]){
        tmp.Sigma = matrix(0, nrow = Ki[v], ncol = Ki[v])
        for (f in 1:(Fl[v] + F0)){
          tmp.Sigma = tmp.Sigma + W_est[[v]][p0, f]*W_est[[v]][p1, f]*tmp.C[[f]]
        }
        
        if (p0==p1){
          diag(tmp.Sigma) = diag(tmp.Sigma) + sigma2_p[[v]][p0]
        }
        
        tmp.Sigma11[(1+(p0-1)*Ki[v]):(p0*Ki[v]), (1+(p1-1)*Ki[v]):(p1*Ki[v])] = tmp.Sigma
      }
    }
    
    Sigma11[k0[v]:k1[v], k0[v]:k1[v]] = tmp.Sigma11
    
  }
  
  ## For off-diagonal elements: v0 != v1 
  for (v0 in 1:V){
    for (v1 in 1:V){
      if (v0!=v1){
        tmp.Sigma11 = matrix(NA, nrow = Ki[v0]*P[v0], ncol = Ki[v1]*P[v1])
        
        tmp.C = vector('list', F0)
        for (f in 1:F0){
          phi.curr0 = matrix(phi_obs[[v0]][[f + Fl[v0]]], nrow = Ki[v0])
          phi.curr1 = matrix(phi_obs[[v1]][[f + Fl[v1]]], nrow = Ki[v1])
          tmp.C[[f]] = phi.curr0 %*% D[[v0]][[f + Fl[v0]]] %*% t(phi.curr1)
        }
        
        for (p0 in 1:P[v0]){
          for (p1 in 1:P[v1]){
            tmp.Sigma = matrix(0, nrow = Ki[v0], ncol = Ki[v1])
            for (f in 1:F0){
              tmp.Sigma = tmp.Sigma + Lambda.est[[v0]][p0, f]*Lambda.est[[v1]][p1, f]*tmp.C[[f]]
            }
            
            tmp.Sigma11[(1+(p0-1)*Ki[v0]):(p0*Ki[v0]), (1+(p1-1)*Ki[v1]):(p1*Ki[v1])] = tmp.Sigma
          }
        }
        
        Sigma11[k0[v0]:k1[v0], k0[v1]:k1[v1]] = tmp.Sigma11
        
      }
    }
  }
  
  # Sigma11 = forceSymmetric(Sigma11)
  
  ## Covariance of Y and zeta ##
  Sigma12_zeta = matrix(NA, nrow = sum(Ki*P), ncol = L0)
  for (v in 1:V){
    
    tmp.Sigma12_zeta = matrix(NA, nrow = Ki[v]*P[v], ncol = L0)
    for (f in (Fl[v] + 1):(Fl[v] + F0)){
      phi.curr = matrix(phi_obs[[v]][[f]], nrow = Ki[v])
      d.curr = diag(D[[v]][[f]])
      
      for (l in 1:L[[v]][f]){
        index.row = 0
        for (p in 1:P[v]){
          for (k in 1:Ki[v]){
            index.row = index.row + 1
            tmp.Sigma12_zeta[index.row, l] = d.curr[l]*phi.curr[k, l]*W_est[[v]][p, f]
          } 
        }  ## End loop of p
        
      } ## End loop of l
    } ## End loop of f 
    
    Sigma12_zeta[k0[v]:k1[v], ] = tmp.Sigma12_zeta
  }
  
  # Covariance of Y and U ##
  Sigma12_U = vector('list', V)
  for (v in 1:V) Sigma12_U[[v]] = matrix(NA, nrow = sum(Ki*P), ncol = Ki[v]*F0)
  
  for (v0 in 1:V){
    for (v1 in 1:V){
      tmp.Sigma12_U = matrix(NA, nrow = Ki[v1]*P[v1], ncol = Ki[v0])
      
      tmp.C = vector('list', F0)
      for (f in 1:F0){
        phi.curr0 = matrix(phi_obs[[v0]][[f + Fl[v0]]], nrow = Ki[v0])
        phi.curr1 = matrix(phi_obs[[v1]][[f + Fl[v1]]], nrow = Ki[v1])
        tmp.C[[f]] = phi.curr1 %*% D[[v1]][[f + Fl[v1]]] %*% t(phi.curr0)
      }
      
      for (f in 1:F0){
        for (p in 1:P[v1]){
          tmp.Sigma12_U[(1+(p-1)*Ki[v1]):(p*Ki[v1]), 
                        (1+(f-1)*Ki[v0]):(f*Ki[v0])] = Lambda.est[[v1]][p, f]*tmp.C[[f]]
        }
      }
      
      Sigma12_U[[v0]][k0[v1]:k1[v1], ] = tmp.Sigma12_U
      
    }
  }
  
  tmp.inv = solve(Sigma11, Yi_all)
  
  zeta_est[i, ] = as.vector(t(Sigma12_zeta) %*% tmp.inv)
  for (v in 1:V){
    U_est[[v]][[i]] = as.vector(t(Sigma12_U[[v]]) %*% tmp.inv)
  }
  
  if (i%%10==0) cat(sprintf("i = %d\n", i))
}  ## End loop of i

#### Conditional expectation to estimate FPC scores and Zi(t) ####
xi_est = vector('list', V)
Z_est = vector('list', V)

for (v in 1:V){
  L1 = L[[v]][1:Fl[v]]
  xi_est[[v]] = matrix(NA, nrow = N, ncol = sum(L1))
  Z_est[[v]] = matrix(NA, nrow = K[v], ncol = Fl[v])
  
  ##
  metadata = metadata_all[[v]]
  Y = Y_all[[v]]
  P = ncol(Y)
  
  ID = ID_all[[v]]
  nobs = table(ID)
  k0 = unname(c(0, cumsum(nobs)[1:(N-1)]) + 1)
  k1 = unname(cumsum(nobs))
  
  D = vector('list', Fl[v] + F0)
  phi_obs = vector('list', Fl[v] + F0)
  for (f in 1:Fl[v]){
    vectors = cov_est_result[[v]]$phi_est[[f]]
    values = cov_est_result[[v]]$d_est[[f]]
    D[[f]] = matrix(0, nrow = L[[v]][f], ncol = L[[v]][f])
    diag(D[[f]]) = values
    
    ## Smooth the eigenfunction on the time of each subject
    phi_obs[[f]] = matrix(NA, nrow = K[v], ncol = L[[v]][f])
    for (l in 1:L[[v]][f]){
      phi_obs[[f]][, l] = bs.smooth(vectors[, l], t.argvals, metadata$time)$est.value
    }
  }
  
  for (f in (Fl[v] + 1):(Fl[v] + F0)){
    vectors = psi_est[[f - Fl[v]]]
    values = d0_est[[f - Fl[v]]]
    D[[f]] = matrix(0, nrow = L[[v]][f], ncol = L[[v]][f])
    diag(D[[f]]) = values
    
    ## Smooth the eigenfunction on the time of each subject
    phi_obs[[f]] = matrix(NA, nrow = K[v], ncol = L[[v]][f])
    for (l in 1:L[[v]][f]){
      phi_obs[[f]][, l] = bs.smooth(vectors[, l], t.argvals, metadata$time)$est.value
    }
  }
  
  W_est = cbind(cov_est_result[[v]]$W_est, Lambda.est[[v]])
  sigma2_p = sigma_p_all[[v]]^2
  
  ## For each subject
  for (i in 1:N){
    Yi = matrix(Y[metadata$ID==i, ], ncol = P)
    Ki = nrow(Yi)
    Yi = as.vector(Yi)  ## vectorize Y: 1st column; 2nd column, etc.
    Yi = matrix(Yi, ncol = 1)
    
    tmp.C = vector('list', Fl[v] + F0)
    for (f in 1:(Fl[v] + F0)){
      phi.curr = matrix(phi_obs[[f]][k0[i]:k1[i], ], nrow = Ki)
      tmp.C[[f]] = phi.curr %*% D[[f]] %*% t(phi.curr)
    }
    
    ## Covariance of Y
    Sigma11 = matrix(NA, nrow = Ki*P, ncol = Ki*P)
    for (p0 in 1:P){
      for (p1 in 1:P){
        tmp.Sigma = matrix(0, nrow = Ki, ncol = Ki)
        for (f in 1:(Fl[v] + F0)){
          tmp.Sigma = tmp.Sigma + W_est[p0, f]*W_est[p1, f]*tmp.C[[f]]
        }
        
        if (p0==p1){
          diag(tmp.Sigma) = diag(tmp.Sigma) + sigma2_p[p0]
        }
        
        Sigma11[(1+(p0-1)*Ki):(p0*Ki), (1+(p1-1)*Ki):(p1*Ki)] = tmp.Sigma
      }
    }
    
    Sigma11 = forceSymmetric(Sigma11)
    
    ## Covariance of Y and xi ##
    Sigma12 = matrix(NA, nrow = Ki*P, ncol = sum(L1))
    index = 0
    for (f in 1:Fl[v]){
      phi.curr = matrix(phi_obs[[f]][k0[i]:k1[i], ], nrow = Ki)
      d.curr = diag(D[[f]])
      for (l in 1:L[[v]][f]){
        index = index + 1
        index.row = 0
        for (p in 1:P){
          for (k in 1:Ki){
            index.row = index.row + 1
            Sigma12[index.row, index] = d.curr[l]*phi.curr[k, l]*W_est[p, f]
          }
        } 
        
      }  ## End loop of l
      
    } ## End loop of f
    
    # Covariance of Y and Z ##
    Sigma12_Z = matrix(NA, nrow = Ki*P, ncol = Ki*Fl[v])
    for (f in 1:Fl[v]){
      phi.curr = matrix(phi_obs[[f]][k0[i]:k1[i], ], nrow = Ki)
      tmp.C = phi.curr %*% D[[f]] %*% t(phi.curr)
      for (p in 1:P){
        Sigma12_Z[(1+(p-1)*Ki):(p*Ki), (1+(f-1)*Ki):(f*Ki)] = W_est[p, f]*tmp.C
      }
    }
    
    tmp.inv = solve(Sigma11, Yi)
    Zi = t(Sigma12_Z) %*% tmp.inv
    for (f in 1:Fl[v]){
      Z_est[[v]][k0[i]:k1[i], f] = Zi[(1+(f-1)*Ki):(f*Ki)]
    }
    
    xi_est[[v]][i, ] = as.vector(t(Sigma12) %*% tmp.inv)
    if (i%%10==0) cat(sprintf("i = %d\n", i))
  }  ## End loop of i
  
}

save(list = c("xi_est", 'zeta_est', 'Z_est', 'U_est'), 
     file = 'RData/new/FPCs_new1.RData')

