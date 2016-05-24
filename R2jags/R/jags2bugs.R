jags2bugs <- function(path=getwd(), parameters.to.save, n.chains=3, n.iter=2000, n.burnin=1000, n.thin=2, DIC=TRUE){
  setwd(path)
  #require(R2WinBUGS)
  fit <- jags.sims(parameters.to.save, n.chains=n.chains, n.iter=n.iter, n.burnin=n.burnin, n.thin=n.thin, DIC = DIC)
  class(fit) <- "bugs"
  return(fit)
}                       
