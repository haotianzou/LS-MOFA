FPCA.long = function(ydata, tnew = NULL, npc = NULL, pve = 0.95){
  uID = unique(ydata[, 1])
  N = length(uID)
  lt = length(tnew)
  newdata = data.frame(t0 = rep(tnew, each = lt), t1 = rep(tnew, lt))
  
  # 
  Yprod = matrix(NA, nrow = 0, ncol = 3)
  colnames(Yprod) = c('t0', 't1', 'Y')
  for (i in 1:N){
    tmp.Y = ydata[which(ydata[, 1]==uID[i]), ]
    tmp.Y = tmp.Y[!duplicated(tmp.Y[, 2]), ]
    
    tmp.time = tmp.Y[, 2]
    tmp.lt = length(tmp.time)
    t0 = rep(tmp.time, each = tmp.lt)
    t1 = rep(tmp.time, tmp.lt)
    
    tmp.data = matrix(NA, nrow = tmp.lt*tmp.lt, ncol = 3)
    colnames(tmp.data) = c('t0', 't1', 'Y')
    tmp.data[, 1] = t0
    tmp.data[, 2] = t1
    
    tmp.data[, 3] = rep(tmp.Y[, 3], each = tmp.lt)*rep(tmp.Y[, 3], tmp.lt)
    
    Yprod = rbind(Yprod, tmp.data)
  }
  
  Yprod = as.data.frame(Yprod)
  Yprod_offdiag = Yprod[which(Yprod[, 1]!=Yprod[, 2]), ] ## extract off diagonal elements
  Yprod_diag = Yprod[which(Yprod[, 1]==Yprod[, 2]), ]
  
  #
  cov_auto = array(NA, dim = c(lt, lt))
  tmp.Yprod_offdiag = Yprod_offdiag
  
  gam.obj = gam(Y ~ te(t0, t1, bs = 'ts'), data = tmp.Yprod_offdiag)
  pred = predict(gam.obj, newdata)
  tmp.cov = matrix(pred, nrow = lt, ncol = lt, byrow = T)
  tmp.cov = (tmp.cov + t(tmp.cov))/2
  cov_auto = tmp.cov
  C_mat = cov_auto
  
  newdata2 = data.frame(t0 = Yprod_diag[, 1], t1 = Yprod_diag[, 2])
  pred2 = predict(gam.obj, newdata2)
  sigma2 = pmax(0.01, mean(Yprod_diag[, 3] - pred2))
  
  ## Estimate eigenvalues and eigenfunctions
  eigen.obj = eigen(C_mat)
    
  values = eigen.obj$values[eigen.obj$values>0]
  if (is.null(npc)){
    pve_data = cumsum(values)/sum(values)
    npc = min(which(pve>0.95))
  }
  
  values = values[1:npc]/lt
  vectors = matrix(eigen.obj$vectors[, 1:npc]*sqrt(lt), ncol = npc)
  
  phi_est = matrix(NA, nrow = lt, ncol = npc)
  d_est = values
  for (l in 1:npc){
    if (sign(vectors[1, l])!=1) vectors[, l] = -vectors[, l]
    phi_est[, l] = vectors[, l]
  }
  
  # Estimate scores
  scores_long = matrix(NA, nrow = N, ncol = npc)
  for (i in 1:N){
    tmp.Y = ydata[which(ydata[, 1]==uID[i]), ]
    Yi = matrix(tmp.Y[, 3], ncol = 1)
    time_i = tmp.Y[, 2]
    Ki = length(Yi)
    
    ##
    vectors = phi_est
    values = d_est
    D = matrix(0, nrow = npc, ncol = npc)
    diag(D) = values
    
    ## Smooth the eigenfunction on the time of each subject
    phi_obs = matrix(NA, nrow = Ki, ncol = npc)
    for (l in 1:npc){
      phi_obs[, l] = bs.smooth(vectors[, l], tnew, time_i)$est.value
    }
    
    phi.curr = phi_obs
    tmp.C = phi.curr %*% D %*% t(phi.curr)
    
    ## 
    Sigma11 = matrix(NA, nrow = Ki, ncol = Ki)
    tmp.Sigma = matrix(0, nrow = Ki, ncol = Ki)
    tmp.Sigma = tmp.C
        
    diag(tmp.Sigma) = diag(tmp.Sigma) + sigma2
    Sigma11 = tmp.Sigma
    
    ## 
    Sigma12 = matrix(NA, nrow = Ki, ncol = npc)
    for (l in 1:npc){
      Sigma12[, l] = phi.curr[, l]*values[l] 
    }
    
    ## 
    tmp.inv = solve(Sigma11, Yi)
    scores_long[i, ] = as.vector(t(Sigma12) %*% tmp.inv)
  }
  
  return(scores_long)
}
