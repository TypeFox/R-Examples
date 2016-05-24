EMbssn <- function(ti,alpha,beta,delta,loglik=F,accuracy = 1e-8,show.envelope="FALSE",iter.max = 500)
{
  #Running the algorithm
  out <- algEMbssn(ti,alpha,beta,delta,loglik,accuracy,show.envelope,iter.max)

  #show result
  cat('\n')
  cat('---------------------------------------------------------\n')
  cat('Birnbaum-Saunders model based on Skew-Normal distribution\n')
  cat('---------------------------------------------------------\n')
  cat('\n')
  cat('Observations =',length(ti))
  cat('\n')
  cat('-----------\n')
  cat('Estimates\n')
  cat('-----------\n')
  cat('\n')
  print(round(out$result$table,5))
  cat('------------------------\n')
  cat('Model selection criteria\n')
  cat('------------------------\n')
  cat('\n')
  critFin  <- c(out$result$loglik, out$result$AIC, out$result$BIC, out$result$HQC)
  critFin  <- round(t(as.matrix(critFin)),digits=3)
  dimnames(critFin) <- list(c("Value"),c("Loglik", "AIC", "BIC","HQC"))
  print(critFin)
  cat('-------\n')
  cat('Details\n')
  cat('-------\n')
  cat('\n')
  cat('Iterations =',out$result$iter)
  cat('\n')
  cat("Processing time =",out$result$time,units(out$result$time))
  cat('\n')
  cat("Convergence =",out$result$convergence)
  cat('\n')
  res            <- list(iter = out$result$iter,criterion = out$result$criterion, alpha=out$result$alpha, beta=out$result$beta, lambda=out$result$lambda, SE=out$result$EP,table = out$result$table,loglik=out$result$loglik, AIC=out$result$AIC, BIC=out$result$BIC, HQC=out$result$HQC, time = out$result$time, convergence = out$result$convergence)
  obj.out        <- list(res = res)
  class(obj.out) <-  "bssn"
  return(obj.out)
}


# EMbssn(ti,alpha0,beta0,delta0,loglik=T,show.envelope=TRUE)
