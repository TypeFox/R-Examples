summary.bkpc <-
function(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), n.burnin = 0, ...){
  
  
  
  n.samples <- dim(object$beta)[1]
  if (n.burnin >= n.samples) stop("error: too many burn-in iterations specified")
  if (n.burnin < 0) n.burnin <- 0
  
  output <- paste("\n", "  BKPC settings:", "\n", "\n",
                  "  Iterations = ", object$n.iter, "\n", 
                  "  Burn-in = ", (n.burnin * object$thin), "\n",
                  "  Thinning interval = ",object$thin, "\n", "\n", sep = "")

  betasum <- cbind(apply(object$beta[(n.burnin + 1) : n.samples, ], 2, mean),  apply(object$beta[(n.burnin + 1) : n.samples, ], 2, sd))
  colnames(betasum) <- c("Mean", "SD")
  rownames(betasum) <- rownames(betasum, do.NULL = FALSE)
  for(i in 1 : dim(object$beta)[2])  rownames(betasum)[i] <- paste("beta[", i,"]", sep = "")
  
  tausum <- cbind(apply(object$tau[(n.burnin + 1) : n.samples, ], 2, mean),  apply(object$tau[(n.burnin + 1) : n.samples, ], 2, sd))
  colnames(tausum) <- c("Mean", "SD")
  rownames(tausum) <- rownames(tausum, do.NULL = FALSE)
  for(i in 1 : dim(object$tau)[2])  rownames(tausum)[i] <- paste("tau[", i,"]", sep = "")

  sigmasum <- cbind(apply(object$sigmasq[(n.burnin + 1) : n.samples, 1,  drop = FALSE], 2, mean),  apply(object$sigmasq[(n.burnin + 1) : n.samples, 1,  drop = FALSE], 2, sd))
  colnames(sigmasum) <- c("Mean", "SD")
  rownames(sigmasum) <- paste("sigmasq", sep = "")
  
  betaq <- t(apply(object$beta[(n.burnin + 1) : n.samples, ], 2, quantile, probs = quantiles))
  rownames(betaq) <- rownames(betaq, do.NULL = FALSE)
  for(i in 1 : dim(object$beta)[2])  rownames(betaq)[i] <- paste("beta[", i,"]", sep = "")
  tauq <- t(apply(object$tau[(n.burnin + 1) : n.samples, ], 2, quantile, probs = quantiles))
  rownames(tauq) <- rownames(tauq, do.NULL = FALSE)
  for(i in 1 : dim(object$tau)[2])  rownames(tauq)[i] <- paste("tau[", i,"]", sep = "")
  sigmasqq <- t(apply(object$sigmasq[(n.burnin + 1) : n.samples, 1,  drop = FALSE], 2, quantile, quantiles))
  rownames(sigmasqq) <- paste("sigmasq", sep = "")

  cat(output)
  cat("\n  1. Empirical mean and standard deviation for each variable: \n \n")
  print(rbind(betasum, tausum, sigmasum))
  cat("\n\n  2. Quantiles for each variable: \n \n") 
  print(rbind(betaq, tauq, sigmasqq))
  
}
