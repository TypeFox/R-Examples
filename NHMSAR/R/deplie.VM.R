deplie.VM <-
function(par) { 

  m = par[1] %% (2*pi)  
  if (length(par)>2) {kappa = c(exp(par[2]),par[3:length(par)])}
  else {kappa = c(exp(par[2]))}
  res <- NULL
  res$m <- m
  res$kappa <- kappa
  return(res)

}
