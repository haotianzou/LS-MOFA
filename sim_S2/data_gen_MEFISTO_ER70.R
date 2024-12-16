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
  psi2 = list(psi11 = function(t) sqrt(2)*cos(2*pi*t), 
              psi12 = function(t) sqrt(2)*sin(2*pi*t))
  psi3 = list(psi11 = function(t) sqrt(2)*cos(3*pi*t), 
              psi12 = function(t) sqrt(2)*sin(3*pi*t))
  psi = list(psi1 = psi1, psi2 = psi2, psi3 = psi3)
  
  ## Shared FPC score
  L0 = rep(unlist(Lu), length(psi))
  zeta = matrix(NA, nrow = N, ncol = sum(L0))
  d0 = list(d01 = c(0.5, 0.2), d02 = c(0.25, 0.1), d03 = c(0.125, 0.05))
  index = 0
  for (f in 1:length(psi)){
    for (l in 1:L0[f]){
      index = index + 1
      tmp.z = rnorm(N, 0, 1)
      tmp.z = unname(scale(tmp.z))*sqrt(d0[[f]][l])
      zeta[, index] = tmp.z
    }
  }
  
  ## Define Lambda 
  Lambda = vector('list', V)
  for (v in 1:V){
    Lambda[[v]] = matrix(NA, nrow = P[v], ncol = length(psi))
    if (v==1) {
      for (f in 1:length(psi)){
        Lambda[[v]][, f] = rnorm(P[v], 0, 0.5)*rbinom(P[v], 1, 0.5)
      }
      Lambda[[v]][1, ] = 1
    } else {
      for (f in 1:length(psi)){
        Lambda[[v]][, f] = rnorm(P[v], 0, 0.5)*rbinom(P[v], 1, 0.5)
      }
    }
  }
  
  Y_all = vector('list', V)
  metadata_all = vector('list', V)
  xi_all = vector('list', V)
  Z_all = W_all = vector('list', V)
  
  for (v in 1:V){
    
    if (Fl[v]==1){
      s0 = 1; s1 = L[[v]]
    } else {
      s0 = c(1, cumsum(L[[v]])[1:(Fl[v] + F0 - 1)] + 1)
      s1 = cumsum(L[[v]])
    }
    
    phi_v = psi
    d_v = d0
    
    Z = matrix(NA, nrow = K[v], ncol = Fl[v] + F0)
    phi_obs = vector('list', Fl[v] + F0) 
    index = 0
    
    ## Shared FPC scores ##
    xi = zeta
    
    ## Compute the eigenfunction 
    for (f in 1:(Fl[v] + F0)){
      phi_obs[[f]] = matrix(NA, nrow = K[v], ncol = L[[v]][f])
      for (l in 1:L[[v]][f]){
        phi_obs[[f]][, l] = phi_v[[f]][[l]](time_all[[v]])
      }
    }
    
    xi_all[[v]] = xi
    
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
    W = Lambda[[v]]
    
    W2 = Varimax(W)$loadings
    W2 = abs(W2)
    
    W = W2
    if (v==1) W[1, ] = 1
    
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
  gamma = rep(0.4, 6)
  logh0 = 0.8
  
  xi = xi_all[[1]]
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

