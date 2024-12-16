data_gen = function(N, V, P, Fl, F0, L, Lu, seed){
  set.seed(seed)
  
  t.argvals = seq(0, 1, by = 0.01)
  lt = length(t.argvals)
  
  ## generate observed time
  time_all = vector('list', V)
  ID_all = vector('list', V)
  K = rep(NA, V)
  for (v in 1:V){
    time = NULL
    ID = NULL
    for (i in 1:N){
      t0 = 0
      t1 = sample(1:10, size = 1)/100
      t2 = t1 + c(1:9)/10
      t3 = c(t0, t1, t2)
      time = c(time, t3)
      ID = c(ID, rep(i, length(t3)))
    }
    time_all[[v]] = time
    ID_all[[v]] = ID
    K[v] = length(time)
  }
  
  psi1 = list(psi11 = function(t) sqrt(2)*cos(pi*t), 
              psi12 = function(t) sqrt(2)*sin(pi*t))
  
  ## Shared FPC score
  L0 = unlist(Lu)
  zeta = matrix(NA, nrow = N, ncol = sum(L0))
  d0 = list(d01 = c(0.5, 0.2))
  index = 0
  for (f in 1:F0){
    for (l in 1:L0){
      index = index + 1
      tmp.z = rnorm(N, 0, 1)
      tmp.z = unname(scale(tmp.z))*sqrt(d0[[f]][l])
      zeta[, index] = tmp.z
    }
  }
  
  ## Define Lambda 
  Lambda = vector('list', V)
  for (v in 1:V){
    Lambda[[v]] = matrix(NA, nrow = P[v], ncol = F0)
    if (v==1) {
      for (f in 1:F0){
        Lambda[[v]][, f] = rnorm(P[v], 0, 0.5)*rbinom(P[v], 1, 0.5)
      }
      Lambda[[v]][1, ] = 1
    } else {
      for (f in 1:F0){
        Lambda[[v]][, f] = rnorm(P[v], 0, 0.5)*rbinom(P[v], 1, 0.5)
      }
    }
  }
  
  phi1 = list(phi11 = function(t) sqrt(2)*sin(pi*t), 
              phi12 = function(t) sqrt(2)*cos(pi*t))
  phi2 = list(phi11 = function(t) sqrt(2)*sin(2*pi*t), 
              phi12 = function(t) sqrt(2)*cos(2*pi*t))
  phi_v1 = list(phi1 = phi1, phi2 = phi2, psi1 = psi1)
  
  phi1 = list(phi11 = function(t) sqrt(2)*sin(pi*t), 
              phi12 = function(t) sqrt(2)*cos(pi*t))
  phi2 = list(phi11 = function(t) sqrt(2)*sin(2*pi*t), 
              phi12 = function(t) sqrt(2)*cos(2*pi*t))
  phi_v2 = list(phi1 = phi1, phi2 = phi2, psi1 = psi1)
  
  phi = list(phi_v1 = phi_v1, phi_v2 = phi_v2)
  
  d = list(d1 = c(3, 0.5), d2 = c(1.5, 0.5), d01 = d0$d01)
  
  Y_all = vector('list', V)
  metadata_all = vector('list', V)
  xi_all = vector('list', V)
  Z_all = W_all = vector('list', V)
  C = vector('list', V)
  C_view = vector('list', V)  ## view-wise covariance: C^(1, v') (s, t), v' = 2, ..., V
  for (v in 1:V){
    
    if (Fl[v]==1){
      s0 = 1; s1 = L[[v]]
    } else {
      s0 = c(1, cumsum(L[[v]])[1:(Fl[v] + F0 - 1)] + 1)
      s1 = cumsum(L[[v]])
    }
    
    phi_v = phi[[v]]
    d_v = d
    
    Z = matrix(NA, nrow = K[v], ncol = Fl[v] + F0)
    xi = matrix(NA, nrow = N, ncol = sum(L[[v]]))
    phi_obs = vector('list', Fl[v] + F0) 
    index = 0
    
    ## Generate view-specific FPC scores ##
    for (f in 1:Fl[v]){
      for (l in 1:L[[v]][f]){
        index = index + 1
        tmp.z = rnorm(N, 0, 1)
        tmp.z = unname(scale(tmp.z))*sqrt(d_v[[f]][l])
        xi[, index] = tmp.z
      }
    }
    
    ## Shared FPC scores ##
    xi[, (index+1):ncol(xi)] = zeta
    
    ## Compute the eigenfunction and covariance matrix
    C1 = vector('list', Fl[v] + F0)  ## true covariance matrix
    if (v>1) C_view[[v]] = vector('list', F0)
    for (f in 1:(Fl[v] + F0)){
      phi_obs[[f]] = matrix(NA, nrow = K[v], ncol = L[[v]][f])
      C1[[f]] = matrix(NA, nrow = lt, ncol = lt)
      phi_true = matrix(NA, nrow = lt, ncol = L[[v]][f])
      for (l in 1:L[[v]][f]){
        phi_obs[[f]][, l] = phi_v[[f]][[l]](time_all[[v]])
        phi_true[, l] = phi_v[[f]][[l]](t.argvals)
      }
      d_mat = matrix(0, nrow = L[[v]][f], ncol = L[[v]][f])
      diag(d_mat) = d_v[[f]]
      C1[[f]] = phi_true %*% d_mat %*% t(phi_true)
      
      ##
      if (v>1) {
        if (f>Fl[v]) C_view[[v]][[f - Fl[v]]] = phi_true %*% d_mat %*% t(phi_true)
      }
    }
    
    xi_all[[v]] = xi
    C[[v]] = C1
    
    ## Compute Z matrix 
    for (k in 1:K[v]){
      tmp.xi = xi[ID_all[[v]][k], ]
      tmp.phi = NULL
      for (f in 1:(Fl[v] + F0)){
        tmp.phi = c(tmp.phi, phi_obs[[f]][k, ])
      }
      
      for (f in 1:(Fl[v] + F0)){
        Z[k, f] = sum(tmp.xi[s0[f]:s1[f]]*tmp.phi[s0[f]:s1[f]])
      }
    }
    
    ## Generate W matrix
    W = matrix(NA, nrow = P[v], ncol = Fl[v] + F0)
    W[, (Fl[v] + 1):(Fl[v] + F0)] = Lambda[[v]]
    for (p in 1:P[v]){
      w = rnorm(Fl[v], 0, 1)
      s = rbinom(Fl[v], 1, 0.5)
      W[p, 1:Fl[v]] = w*s
    }
    
    W2 = Varimax(W[, 1:Fl[v]])$loadings
    W2 = abs(W2)
    
    W = cbind(W2, Lambda[[v]])
    
    tmp.Y = Z %*% t(W)
    Z_all[[v]] = Z
    W_all[[v]] = W
    
    sigma_p = runif(P[v], 0.6, 1.0)
    err = matrix(NA, K[v], P[v])
    for (p in 1:P[v]){
      err[, p] = rnorm(K[v], 0, sigma_p[p])
    }
    
    Y = tmp.Y + err
    Y_all[[v]] = Y
    
    metadata = data.frame(ID = ID_all[[v]], time = time_all[[v]])
    metadata_all[[v]] = metadata
    
  }
  
  
  ## Generate the survival outcome
  gamma = rep(0.5, 10)
  logh0 = 0
  
  xi = cbind(xi_all[[1]][, 1:4], xi_all[[2]][, 1:4], zeta)
  C_time = runif(N, 0, 3)
  C_time = pmin(C_time, rep(1, N))
  
  surv.dat = data.frame(surv_time = NA, status = NA, xi = xi)
  
  for (i in 1:N){
    Si <- runif(1, 0, 1) ## survival probability
    xi_i = xi[i, ]
    hi <- exp(logh0 + sum(gamma*xi_i))
    ti <- -log(Si)/hi
    
    surv.dat$surv_time[i] <- min(ti, C_time[i]) 
    surv.dat$status[i] <- as.integer(ti<C_time[i])
  }
  
  
  truncate = vector('list', V)
  for (v in 1:V){
    tmp.surv_time = rep(surv.dat$surv_time, each = 11)
    truncate[[v]] = as.numeric(metadata_all[[v]]$time<=tmp.surv_time)
  }
  
  #### Truncate the outcome
  for (v in 1:V){
    Y_all[[v]] = Y_all[[v]][truncate[[v]]==1, ]
    metadata_all[[v]] = metadata_all[[v]][truncate[[v]]==1, ]
    Z_all[[v]] = Z_all[[v]][truncate[[v]]==1, ]
  }
  
  
  l = list(Y = Y_all, metadata = metadata_all, xi = xi_all, W = W_all, Z = Z_all, surv.dat = surv.dat)
  return(l)
}

