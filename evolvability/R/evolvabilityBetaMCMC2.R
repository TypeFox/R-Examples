#' @export
evolvabilityBetaMCMC2 = function(G_mcmc, Beta_mcmc, post.dist=FALSE){
  G_Beta = cbind(G_mcmc, Beta_mcmc)
  dimG = sqrt(ncol(G_mcmc))
  X1 = t(apply(G_Beta, 1,
               function(GB){
                 G = matrix(GB[1:(dimG^2)], ncol=dimG)
                 B = cbind(GB[(dimG^2+1):ncol(G_Beta)])
                 eB = t(B)%*%G%*%B 
                 rB = sqrt(t(B)%*%(G%*%G)%*%B)
                 cB = 1/(t(B)%*%solve(G)%*%B)
                 aB = cB/eB
                 iB = 1-aB
                 c(eB = eB, rB = rB, cB = cB, aB = aB, iB = iB)
               }))
  X = list()
  X$Beta.median = cbind(median = apply(Beta_mcmc, 2, median), coda::HPDinterval(coda::mcmc(Beta_mcmc)))
  X$summary = cbind(median = apply(X1, 2, median), coda::HPDinterval(coda::mcmc(X1)))
  if(post.dist == TRUE){
    X$post.dist = X1
  }
  X$call = match.call()
  class(X) = "evolvabilityBetaMCMC2"
  X
}

#' @export
print.evolvabilityBetaMCMC2 = function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nEvolvability parameters, posterior medians and 95% HPD intervals:\n")
  print(x$summary)
} 
