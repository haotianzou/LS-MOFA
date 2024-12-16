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
library(survival)
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
Fl[v] = 2 + F0  ## 40.3% total variation
# fa.obj = fa(Y_base, nfactors = Fl[v] + F0, rotate = 'varimax')
# loadings = matrix(fa.obj$loadings, ncol = Fl[v] + F0)

## Second omics
v = 2
Y_base = Y_all[[v]][metadata_all[[v]]$time==0, ]
#scree(Y_base)
Fl[v] = 6 + F0  # 44.5% total variation
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

#### Compute W and C  ####
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
  
  ## Baseline information
  Y_base = Y[metadata$time==0, ]
  cov.mat = cov(Y_base)
  
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
  sigma_p_all = sqrt(diag(Psi))
  
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
      tmp.cov1 = cov_auto[, t0, t1] 
      tmp.cov2 = cov_cross[2:P, t0, t1] 
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
  
  cov_est_result = list(phi_est = phi_est, d_est = d_est, W_est = W.curr, L = L)
}

## Conditional expectation to estimate FPC scores and/or Zi(t) ##
L = vector('list', V)
xi_est = vector('list', V)
Z_est = vector('list', V)

for (v in 1:1){
  L1 = cov_est_result$L
  L = L1
  xi_est[[v]] = matrix(NA, nrow = N, ncol = sum(L1))
  Z_est[[v]] = matrix(NA, nrow = K[v], ncol = Fl[v])
  
  metadata = metadata_all[[v]]
  Y = Y_all[[v]]
  P = ncol(Y)
  
  ID = ID_all[[v]]
  nobs = table(ID)
  k0 = unname(c(0, cumsum(nobs)[1:(N-1)]) + 1)
  k1 = unname(cumsum(nobs))
  
  D = vector('list', Fl[v])
  phi_obs = vector('list', Fl[v])
  for (f in 1:Fl[v]){
    vectors = cov_est_result$phi_est[[f]]
    values = cov_est_result$d_est[[f]]
    D[[f]] = matrix(0, nrow = L[f], ncol = L[f])
    diag(D[[f]]) = values
    
    ## Smooth the eigenfunction on the time of each subject
    phi_obs[[f]] = matrix(NA, nrow = K[v], ncol = L[f])
    for (l in 1:L[f]){
      phi_obs[[f]][, l] = bs.smooth(vectors[, l], t.argvals, metadata$time)$est.value
    }
  }
  
  W_est = cov_est_result$W_est
  sigma2_p = sigma_p_all^2
  
  ## For each subject
  for (i in 1:N){
    Yi = matrix(Y[metadata$ID==i, ], ncol = P)
    Ki = nrow(Yi)
    Yi = as.vector(Yi)  ## vectorize Y: 1st column; 2nd column, etc.
    Yi = matrix(Yi, ncol = 1)
    
    ## Covariance of Y
    Sigma11 = matrix(NA, nrow = Ki*P, ncol = Ki*P)
    for (p0 in 1:P){
      for (p1 in 1:P){
        tmp.Sigma = matrix(0, nrow = Ki, ncol = Ki)
        for (f in 1:Fl[v]){
          phi.curr = matrix(phi_obs[[f]][k0[i]:k1[i], ], nrow = Ki)
          tmp.C = phi.curr %*% D[[f]] %*% t(phi.curr)
          tmp.Sigma = tmp.Sigma + W_est[p0, f]*W_est[p1, f]*tmp.C
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
      for (l in 1:L[f]){
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
    
    Zi = t(Sigma12_Z) %*% solve(Sigma11, Yi)
    for (f in 1:Fl[v]){
      Z_est[[v]][k0[i]:k1[i], f] = Zi[(1+(f-1)*Ki):(f*Ki)]
    }
    
    xi_est[[v]][i, ] = as.vector(t(Sigma12) %*% solve(Sigma11, Yi))
    if (i%%10==0) cat(sprintf("i = %d\n", i))
  }  ## End loop of i
  
} ## End loop of v

xi_est_lipid = xi_est[[1]]
Z_est_lipid = Z_est[[1]]

save(list = c("xi_est_lipid", 'Z_est_lipid'), file = 'RData/new/FPCs_lipidomics_new.RData')

