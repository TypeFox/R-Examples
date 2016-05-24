deplie.c.VM <-
function(par,constr) { 

  m = par[1] %% (2*pi) 
  order = length(par)/2-1 
  if (constr) { kappa = c(exp(par[2]),par[3:(length(par)-1)]*exp(1i*par[length(par)]))} else {
  kappa = c(exp(par[2]),par[3:(3+order-1)]+1i*par[(3+order):length(par)])
  }
  res <- NULL
  res$m <- m
  res$kappa <- kappa
  return(res)

}
