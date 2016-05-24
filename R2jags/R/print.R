print.rjags <- function(x, digits = 3, 
  intervals = c(0.025, 0.25, 0.5, 0.75, 0.975), ...)
{
  x <- x$BUGSoutput
  sims.matrix <- x$sims.matrix
  mu.vect <- apply(sims.matrix, 2, mean)
  sd.vect <- apply(sims.matrix, 2, sd)
  int.matrix <- apply(sims.matrix, 2, quantile, intervals)
  if (x$n.chains>1) {
    n.eff <- x$summary[,"n.eff"]
    Rhat <- x$summary[,"Rhat"] 
  } else {
    n.eff <- Rhat <- NULL
  }
  summaryMatrix <- t(rbind(mu.vect, sd.vect, int.matrix, Rhat, n.eff))

  rownameMatrix <- rownames(summaryMatrix)
  dev.idx <- match("deviance", rownameMatrix)
  if(any(!is.na(dev.idx))){
    summaryMatrix <- rbind(summaryMatrix[-dev.idx,], summaryMatrix[dev.idx,])
    rownames(summaryMatrix) <- c(rownameMatrix[-dev.idx], rownameMatrix[dev.idx])
  }
  
  if (!is.null(x$model.file)) 
      cat("Inference for Bugs model at \"", x$model.file, "\", ", 
          sep = "")
  if (!is.null(x$program)) 
      cat("fit using ", x$program, ",", sep = "")
  cat("\n ", x$n.chains, " chains, each with ", x$n.iter, " iterations (first ", 
      x$n.burnin, " discarded)", sep = "")
  if (x$n.thin > 1) 
      cat(", n.thin =", x$n.thin)
  cat("\n n.sims =", x$n.sims, "iterations saved\n")
  print(round(summaryMatrix, digits), ...)
  if (x$n.chains > 1) {
      cat("\nFor each parameter, n.eff is a crude measure of effective sample size,")
      cat("\nand Rhat is the potential scale reduction factor (at convergence, Rhat=1).\n")
  }
  if (x$isDIC) {
      msgDICRule <- ifelse(x$DICbyR, "(using the rule, pD = var(deviance)/2)", 
          "(using the rule, pD = Dbar-Dhat)")
      cat(paste("\nDIC info ", msgDICRule, "\n", sep = ""))
      if (length(x$DIC) == 1) {
          cat("pD =", fround(x$pD, 1), "and DIC =", fround(x$DIC, 
              1))
      }
      else if (length(x$DIC) > 1) {
          print(round(x$DIC, 1))
      }
      cat("\nDIC is an estimate of expected predictive error (lower deviance is better).\n")
  }
  invisible(x)
  #print(x$BUGSoutput,...)
} 




#function (x, digits.summary = 1, ...) 
#{
#    if (!is.null(x$model.file)) 
#        cat("Inference for Bugs model at \"", x$model.file, "\", ", 
#            sep = "")
#    if (!is.null(x$program)) 
#        cat("fit using ", x$program, ",", sep = "")
#    cat("\n ", x$n.chains, " chains, each with ", x$n.iter, " iterations (first ", 
#        x$n.burnin, " discarded)", sep = "")
#    if (x$n.thin > 1) 
#        cat(", n.thin =", x$n.thin)
#    cat("\n n.sims =", x$n.sims, "iterations saved\n")
#    print(round(x$summary, digits.summary), ...)
#    if (x$n.chains > 1) {
#        cat("\nFor each parameter, n.eff is a crude measure of effective sample size,")
#        cat("\nand Rhat is the potential scale reduction factor (at convergence, Rhat=1).\n")
#    }
#    if (x$isDIC) {
#        msgDICRule <- ifelse(x$DICbyR, "(using the rule, pD = var(deviance)/2)", 
#            "(using the rule, pD = Dbar-Dhat)")
#        cat(paste("\nDIC info ", msgDICRule, "\n", sep = ""))
#        if (length(x$DIC) == 1) {
#            cat("pD =", fround(x$pD, 1), "and DIC =", fround(x$DIC, 
#                1))
#        }
#        else if (length(x$DIC) > 1) {
#            print(round(x$DIC, 1))
#        }
#        cat("\nDIC is an estimate of expected predictive error (lower deviance is better).\n")
#    }
#    invisible(x)
#}
#
