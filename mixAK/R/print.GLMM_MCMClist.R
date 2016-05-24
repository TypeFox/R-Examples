##
##  PURPOSE:   Generalized linear mixed model with possibly several response variables
##             and normal mixtures in the distribution of continuous response and random effects
##             * two parallel chains
##             * print method for the resulting object
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:    04/11/2011
##
##  FUNCTIONS:  print.GLMM_MCMClist
##
## ================================================================================================

## *************************************************************
## print.GLMM_MCMClist
## *************************************************************
print.GLMM_MCMClist <- function(x, ...)
{
  R <- sum(x[[1]]$R)
  chNames <- c("Chain 1", "Chain 2")
  
  cat("\n")
  cat(paste("     Generalized linear mixed model for ", R, " responses estimated using MCMC\n", sep=""))
  cat("     ====================================================================\n")

  cat("\nPenalized expected deviance:")
  cat("\n----------------------------\n")
  print(x$PED, ...)

  #cat("\nDeviance information criteria:")
  #cat("\n------------------------------\n")
  #xDIC <- rbind(x[[1]]$DIC, x[[2]]$DIC)
  #rownames(xDIC) <- chNames
  #print(xDIC, ...)  
  
  cat("\nDeviance posterior summary statistics:")
  cat("\n-----------------------------------------------\n")
  pr.Dev <- rbind(x[[1]]$summ.Deviance[, "Deviance"], x[[2]]$summ.Deviance[, "Deviance"])
  rownames(pr.Dev) <- chNames
  colnames(pr.Dev) <- rownames(x[[1]]$summ.Deviance)  
  print(pr.Dev, ...)
  
  if (x[[1]]$lalpha){
    cat("\nPosterior summary statistics for fixed effects:")
    cat("\n-----------------------------------------------\n")
    if (x[[1]]$lalpha == 1){
      pr.alpha <- rbind(x[[1]]$summ.alpha, x[[2]]$summ.alpha)
      rownames(pr.alpha) <- chNames
      colnames(pr.alpha) <- names(x[[1]]$summ.alpha)  
      print(pr.alpha, ...)      
    }else{
      pr.alpha <- cbind(x[[1]]$summ.alpha[,1], x[[2]]$summ.alpha[,1])
      for (i in 2:x[[1]]$lalpha) pr.alpha <- cbind(pr.alpha, x[[1]]$summ.alpha[,i], x[[2]]$summ.alpha[,i])
      rownames(pr.alpha) <- rownames(x[[1]]$summ.alpha)
      colnames(pr.alpha) <- paste(rep(colnames(x[[1]]$summ.alpha), each=2), "(Chain ", rep(1:2, x[[1]]$lalpha), ")", sep="")
      print(pr.alpha, ...)      
    }  
  }

  if (x[[1]]$dimb){
    if (x[[1]]$prior.b$priorK == "fixed"){
      cat(paste("\nDistribution of random effects is a normal mixture with ", x[[1]]$prior.b$Kmax, " components", sep=""))
      cat("\n---------------------------------------------------------------------")
    }else{
      cat(paste("\nDistribution of random effects is a normal mixture with at most ", x$prior.b$Kmax, " components", sep=""))
      cat("\n------------------------------------------------------------------------------------------------")
      cat("\nPosterior distribution of K (number of mixture components):")
      cat("\n-----------------------------------------------------------")
      Kavail <- unique(c(names(x[[1]]$propK_b), names(x[[2]]$propK_b)))
      xpropK_b<- matrix(0, nrow=2, ncol=length(Kavail))
      rownames(xpropK_b) <- chNames
      colnames(xpropK_b) <- Kavail
      xpropK_b[1, names(x[[1]]$propK_b)] <- x[[1]]$propK_b
      xpropK_b[2, names(x[[2]]$propK_b)] <- x[[2]]$propK_b
      print(xpropK_b, ...)
    }
    cat("\nPosterior summary statistics for moments of mixture for random effects:")
    cat("\n-----------------------------------------------------------------------")
    if (x[[1]]$dimb == 1){
      cat("\nMean:\n")
      pr.b.Mean <- rbind(x[[1]]$summ.b.Mean, x[[2]]$summ.b.Mean)
      rownames(pr.b.Mean) <- chNames
      colnames(pr.b.Mean) <- names(x[[1]]$summ.b.Mean)  
      print(pr.b.Mean, ...)      

      cat("\nStandard deviation:\n")
      pr.b.SDCorr <- rbind(x[[1]]$summ.b.SDCorr, x[[2]]$summ.b.SDCorr)
      rownames(pr.b.SDCorr) <- chNames
      colnames(pr.b.SDCorr) <- names(x[[1]]$summ.b.SDCorr)  
      print(pr.b.SDCorr, ...)      
    }else{
      cat("\nMeans:\n")
      pr.b.Mean <- cbind(x[[1]]$summ.b.Mean[,1], x[[2]]$summ.b.Mean[,1])
      for (i in 2:x[[1]]$dimb) pr.b.Mean <- cbind(pr.b.Mean, x[[1]]$summ.b.Mean[,i], x[[2]]$summ.b.Mean[,i])
      rownames(pr.b.Mean) <- rownames(x[[1]]$summ.b.Mean)
      colnames(pr.b.Mean) <- paste(rep(colnames(x[[1]]$summ.b.Mean), each=2), "(Chain ", rep(1:2, x[[1]]$dimb), ")", sep="")
      print(pr.b.Mean, ...)      

      cat("\nStandard deviations and correlations:\n")
      pr.b.SDCorr <- cbind(x[[1]]$summ.b.SDCorr[,1], x[[2]]$summ.b.SDCorr[,1])
      for (i in 2:ncol(x[[1]]$summ.b.SDCorr)) pr.b.SDCorr <- cbind(pr.b.SDCorr, x[[1]]$summ.b.SDCorr[,i], x[[2]]$summ.b.SDCorr[,i])
      rownames(pr.b.SDCorr) <- rownames(x[[1]]$summ.b.SDCorr)
      colnames(pr.b.SDCorr) <- paste(rep(colnames(x[[1]]$summ.b.SDCorr), each=2), "(Chain ", rep(1:2, ncol(x[[1]]$summ.b.SDCorr)), ")", sep="")
      print(pr.b.SDCorr, ...)      
    }      
  }

  if (x[[1]]$R["Rc"]){
    cat("\nPosterior summary statistics for standard deviations")
    cat("\nof residuals of continuous responses:")
    cat("\n----------------------------------------------------\n")
    if (x[[1]]$R["Rc"] == 1){
      pr.sigma_eps <- rbind(x[[1]]$summ.sigma_eps, x[[2]]$summ.sigma_eps)
      rownames(pr.sigma_eps) <- chNames
      colnames(pr.sigma_eps) <- names(x[[1]]$summ.sigma_eps)  
      print(pr.sigma_eps, ...)      
    }else{
      pr.sigma_eps <- cbind(x[[1]]$summ.sigma_eps[,1], x[[2]]$summ.sigma_eps[,1])
      for (i in 2:x[[1]]$R["Rc"]) pr.sigma_eps <- cbind(pr.sigma_eps, x[[1]]$summ.sigma_eps[,i], x[[2]]$summ.sigma_eps[,i])
      rownames(pr.sigma_eps) <- rownames(x[[1]]$summ.sigma_eps)
      colnames(pr.sigma_eps) <- paste(rep(colnames(x[[1]]$summ.sigma_eps), each=2), "(Chain ", rep(1:2, x[[1]]$R["Rc"]), ")", sep="")
      print(pr.sigma_eps, ...)      
    }  
  }  
    
  return(invisible(x))  
}
