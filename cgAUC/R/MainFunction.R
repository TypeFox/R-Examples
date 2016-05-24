# ------------------------------------------------

optimal.delta = function(y, z, l, h, ind.d.l){
  l.i = matrix(rep(l, times = 50), nrow = 50, byrow = TRUE)
  
  delta = seq(0, 5, length = 50)
  m = delta %*% t(ind.d.l)
  
  l.i = l.i + m
  
  l.i.max = apply(l.i, 1, max)
  l.i = l.i / l.i.max
  
  theta = rep(0, 50) # AUC = theta
  
  for(i in 2:50){
    theta[i] = c_cntin(y, z, l.i[i, ], h)[1]
  }
  
  delta.star = delta[which(theta == max(theta))]
  
  return(delta.star)
}

# ------------------------------------------------

cgAUC = function(x, z, h, delta = 1, auto = FALSE, tau = 1, scale = 1){
  if(scale == 0){
    x = as.matrix(x);
    z = as.matrix(z);
  }
  else{
    x = scale(x);
    z = scale(z);
  }
  
  conv = FALSE
  n = dim(x)[1]
  p = dim(x)[2]
  cntin.ri = dscrt.ri = rep(0, p)
  id = diag(p)
  
  for(i in 1:p){
    dscrt.ri[i] = c_dscrt(x, z, id[i, ]   )[1]
    cntin.ri[i] = c_cntin(x, z, id[i, ], h)[1]
  }
  
  beta.i = ifelse(cntin.ri > 0.5, 1, -1)
  
  dscrt.ri = ifelse(dscrt.ri > 0.5, dscrt.ri, (1 - dscrt.ri))
  cntin.ri = ifelse(cntin.ri > 0.5, cntin.ri, (1 - cntin.ri))
  
  y = x * matrix(beta.i, n, p, byrow = TRUE)
  max.x = which(cntin.ri == max(cntin.ri))
  
  theta.sh.h.p = 0
  l = id[max.x, ]
  
  # TGDM
  while(conv == FALSE){
    # Step 1
    d.l = c_d_theta_sh_h_p(y, z, l, h)
    
    max.d.l = max(abs(d.l))
    ind.d.l = ifelse(abs(d.l) >= (tau * max.d.l), 1, 0) * d.l
    # Step 3
    if (auto == TRUE){
      delta = optimal.delta(y, z, l, h, ind.d.l)
    }
    
    l = l + delta * ind.d.l
    l = l / max(l)
    theta.temp = c_cntin(y, z, l, h)[1]
    
    ifelse(abs(theta.temp - theta.sh.h.p) < 0.0001, conv <- TRUE, conv <- FALSE)
    theta.sh.h.p = theta.temp
  }
  
  optimal.dscrt = c_dscrt(y, z, l)
  theta.sh.h.p.var = c_cntin(y, z, l, h)[2]
  
  l = l * beta.i
  
  if(min(l) == -1){
    l = l * -1
    Rev = 1
  }
  else{
    Rev = 0
  }
  
  return(list(
    Rev = Rev,
    l = l,
    
    theta.sh.h.p = theta.sh.h.p,
    theta.sh.h.p.var = theta.sh.h.p.var,
    cntin.ri = cntin.ri,
    
    theta.h.p = optimal.dscrt[1],
    theta.h.p.var = optimal.dscrt[2],
    dscrt.ri = dscrt.ri,
    
    delta = delta
  ))
}
