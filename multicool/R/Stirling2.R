Stirling2 = function(n, k){
  if(any(n < 0) || any(k < 0))
    stop("n and k must be positive integers, n >= 0, k >= 0")
  
  nN = length(n)
  nK = length(k)
  
  if(nN > 1 || nK > 1){
    if(nN != nK){
      grid = expand.grid(n = n, k = k)
      n = grid$n
      k = grid$k
      
      nN = nK = nrow(grid)
    }
  }
  
  result = rep(0, nN)
  for(r in 1:nN){
    result[r] = Stirling2C(n[r], k[r])
  }
  
  return(result)
}

S2 = function(n, k){Stirling2(n, k)}

Bell = function(n){
  if(any(n <=0))
    stop("n must be greater than or equal to 1")
  
  nN = length(n)
  result = rep(0, nN)
  
  for(r in 1:nN)
    result[r] = sum(Stirling2(n[r], 1:n[r]))
  
  return(result)
}

B = function(n){Bell(n)}
