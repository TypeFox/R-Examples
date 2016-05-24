LS.kalman = function(series, start, order = c(p = 0, q = 0), ar.order = NULL, ma.order = NULL, sd.order = NULL, d.order = NULL, include.d = FALSE, m = NULL){
  
  x = start
  T. = length(series)
  
  if(is.null(ar.order)){ar.order = rep(0,order[1])}
  if(is.null(ma.order)){ma.order = rep(0,order[2])}
  if(is.null(sd.order)){sd.order = 0}
  if(is.null(d.order)){d.order = 0}
  if(is.null(m)){m = trunc(0.25*T.^0.8)}
  
  M = m+1
  u = (1:T.)/T.
  
  p = na.omit(c(ar.order, ma.order, sd.order))
  if(include.d == TRUE){
    p = na.omit(c(ar.order, ma.order, d.order, sd.order))
  }
  
  phi. = numeric()
  theta. = numeric()
  sigma. = numeric()
  d. = numeric()
  
  for(j in 1:length(u)){
    X = numeric()
    k = 1
    for(i in 1:length(p)){
      X[i] = sum(x[k:(k+p[i])]*u[j]^(0:p[i]))
      k = k+p[i]+1
    }
    
    phi = numeric()
    k = 1
    if(order[1]>0){
      phi[is.na(ar.order) == 1] = 0
      phi[is.na(ar.order) == 0] = X[k:(length(na.omit(ar.order)))]	
      k = length(na.omit(ar.order))+1
      phi. = rbind(phi., phi)
    }
    
    
    theta = numeric()
    if(order[2]>0){
      theta[is.na(ma.order) == 1] = 0
      theta[is.na(ma.order) == 0] = X[k:(length(na.omit(ma.order))+k-1)]	
      k = length(na.omit(ma.order))+k
      theta. = rbind(theta., theta)
    }
    
    d = 0
    if(include.d == TRUE){
      d = X[k]		
      k = k + 1
      d. = c(d.,d)
    }
    
    sigma = X[k]	
    sigma. = c(sigma., sigma)
  }
  
  sigma = sigma.
  
  Omega = matrix(0, nrow = M, ncol = M)
  diag(Omega) = 1
  
  X = rep(0, M)
  delta = vector("numeric") 
  hat.y = vector("numeric") 
  for(i in 1:T.){
    if(is.null(dim(phi.))==1 & is.null(dim(theta.))==1){psi = c(1,ARMAtoMA(ar = numeric(), ma = numeric(), lag.max = m))}
    if(is.null(dim(phi.))==1 & is.null(dim(theta.))==0){psi = c(1,ARMAtoMA(ar = numeric(), ma = theta.[i,], lag.max = m))}
    if(is.null(dim(phi.))==0 & is.null(dim(theta.))==1){psi = c(1,ARMAtoMA(ar = phi.[i,], ma = numeric(), lag.max = m))}
    if(is.null(dim(phi.))==0 & is.null(dim(theta.))==0){psi = c(1,ARMAtoMA(ar = phi.[i,], ma = theta.[i,], lag.max = m))}
    psi. = numeric()
    if(include.d == TRUE){
      eta = gamma(0:m+d.[i])/(gamma(0:m+1)*gamma(d.[i]))
      for(k in 0:m){
        psi.[k+1] = sum(psi[1:(k+1)]*rev(eta[1:(k+1)]))
      }
      psi = psi.
    }
    g = sigma[i]*rev(psi)
    aux = Omega %*% g
    delta[i] = g %*% aux
    F. = matrix(0,M-1,M-1)
    diag(F.)=1
    F.=cbind(0,F.)
    F.=rbind(F.,0)
    Theta = c(F.%*%aux)
    Q = matrix(0,M,M)
    Q[M,M] = 1
    if(is.na(series[i])){
      Omega = F.%*%Omega%*%t(F.) + Q
      hat.y[i] = t(g)%*%X
      X = F.%*%X
    }
    else{
      Omega = F.%*%Omega%*%t(F.) + Q - Theta %*% solve(delta[i]) %*% Theta
      hat.y[i] = t(g)%*%X
      X = F.%*%X + Theta %*% solve(delta[i])%*%(series[i]-hat.y[i])
    }
  }
  residuals = (series-hat.y)/sqrt(delta[1:T.])
  fitted.values = hat.y
  return(list(residuals = residuals, fitted.values = fitted.values, delta = delta))
}