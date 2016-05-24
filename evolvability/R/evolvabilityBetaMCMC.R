#' @export

evolvabilityBetaMCMC = function(G_mcmc, Beta, post.dist = FALSE){
  X1 = t(apply(G_mcmc, 1,
               function(G){
                 G = matrix(G, ncol=sqrt(length(G)))
                   eB = diag(t(Beta)%*%G%*%Beta) 
                   rB = sqrt(diag(t(Beta)%*%(G%*%G)%*%Beta))
                   cB = 1/diag(t(Beta)%*%solve(G)%*%Beta)
                   aB = cB/eB
                   iB = 1-aB
                 c(eB = eB, rB = rB, cB = cB, aB = aB, iB = iB)
               }))
  n = ncol(Beta)
  X2 = list(eB = X1[,1:n], rB = X1[,(n+1):(2*n)], cB = X1[,(2*n+1):(3*n)], 
            aB = X1[,(3*n+1):(4*n)], iB = X1[,(4*n+1):(5*n)])
  X_summary = cbind(sapply(X2, function(x) apply(x, 1, mean)))
  colnames(X_summary) = c("e_mean", "r_mean", "c_mean", "a_mean", "i_mean")
  X = lapply(X2, function(x) rbind(median = apply(x, 2, median), t(coda::HPDinterval(coda::mcmc(x)))))
  X$Beta = Beta
  X$summary = rbind(median = apply(X_summary, 2, median), t(coda::HPDinterval(coda::mcmc(X_summary))))
  X$call = match.call()
  class(X) = "evolvabilityBetaMCMC"
  X2$summary = X_summary
  if(post.dist) X$post.dist = X2
  X
}

#' @export
print.evolvabilityBetaMCMC = function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nEvolvability, posterior medians and 95% HPD intervals:\n")
  print(t(x$eB))
  cat("\nRespondability, posterior medians and 95% HPD intervals:\n")
  print(t(x$rB))
  cat("\nConditional evolvability, posterior medians and 95% HPD intervals:\n")
  print(t(x$cB))
  cat("\nAutonomy, posterior medians and 95% HPD intervals:\n")
  print(t(x$aB))
  cat("\nIntegration, posterior medians and 95% HPD intervals:\n")
  print(t(x$iB))
} 

#' @export
summary.evolvabilityBetaMCMC = function(object, ...){
  X = list()
  X$call = object$call
  X$Averages = object$summary
  X$Minimum = sapply(object[1:5], function(x) x[, which(x[1,]==min(x[1,]))])
  colnames(X$Minimum) = paste(colnames(X$Min), "_min", sep="")
  X$Maximum = sapply(object[1:5], function(x) x[, which(x[1,]==min(x[1,]))])
  colnames(X$Maximum) = paste(colnames(X$Max), "_max", sep="")
  class(X) = "summary.evolavbilityBetaMCMC"
  X
}

#' @export
print.summary.evolvabilityBetaMCMC = function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nAverage:\n")
  print(x$Averages)
  cat("\nMinimum (the direction with the lowest posterior median):\n")
  print(x$Minimum)
  cat("\nMaximum (the direction with the highest posterior median):\n")
  print(x$Maximum)
}
