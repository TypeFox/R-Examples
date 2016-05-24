LS.whittle.loglik = function(x, series, order = c(p = 0, q = 0), ar.order = NULL, ma.order = NULL, sd.order = NULL, d.order = NULL, include.d = FALSE, N = NULL, S = NULL, include.taper = TRUE){
  
  y = series
  T. = length(y)
  
  if(is.null(N)){N = trunc(T.^0.8)}
  if(is.null(S)){S = trunc(0.2*N)}
  
  M = trunc((T. - N)/S + 1)
  
  if(is.null(ar.order)){ar.order = rep(0,order[1])}
  if(is.null(ma.order)){ma.order = rep(0,order[2])}
  if(is.null(sd.order)){sd.order = 0}
  if(is.null(d.order)){d.order = 0}
  
  
  p = na.omit(c(ar.order, ma.order, sd.order))
  if(include.d == TRUE){
    p = na.omit(c(ar.order, ma.order, d.order, sd.order))
  }
  
  if(length(x) != sum(p+1)){
    stop("error in the number of parameters")	
  }
  
  else{
    
    lik = 0
    for(j in 1:M){
      u = (N/2 + S*(j-1))/T.
      aux = periodogram(y[(1 + S*(j-1)):(N + S*(j-1))], include.taper = TRUE, plot = FALSE)
      I = aux$periodogram
      
      X = numeric()
      k = 1
      for(i in 1:length(p)){
        X[i] = sum(x[k:(k+p[i])]*u^(0:p[i]))
        k = k+p[i]+1
      }
      
      phi = numeric()
      k = 1
      if(order[1]>0){
        phi[is.na(ar.order) == 1] = 0
        phi[is.na(ar.order) == 0] = X[k:(length(na.omit(ar.order)))]	
        k = length(na.omit(ar.order))+1
      }
      
      theta = numeric()
      if(order[2]>0){
        theta[is.na(ma.order) == 1] = 0
        theta[is.na(ma.order) == 0] = X[k:(length(na.omit(ma.order))+k-1)]	
        k = length(na.omit(ma.order))+k
      }
      
      d = 0
      if(include.d == TRUE){
        d = X[k]		
        k = k + 1
      }
      
      sigma = X[k]	
      
      f = fdensity(ar = phi, ma = theta, d = d, sd = sigma, lambda = aux$lambda)
      
      lik = sum(log(f) + I/f)/N + lik
    }
    
    lik = lik/M
    
    if(is.na(sigma) | sigma <= 0){lik = Inf}
    
    lik
    
  }
}