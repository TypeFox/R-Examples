##
##  PURPOSE:   Reversible jump MCMC for a normal mixture model
##             * two parallel chains
##             * print method for the resulting object
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   12/11/2008
##
##  FUNCTION:  print.NMixMCMClist 
##             
##
## ======================================================================
## *************************************************************
## print.NMixMCMC
## *************************************************************
print.NMixMCMClist <- function(x, ped, dic, ...)
{
  cat("\n")
  chNames <- c("Chain 1", "Chain 2")

  if (missing(dic)) dic <- (x[[1]]$prior$priorK == "fixed")
  if (!is.logical(dic)) stop("dic must be logical")

  if (missing(ped)) ped <- (x[[1]]$prior$priorK == "fixed")
  if (!is.logical(ped)) stop("dic must be logical")
  
  if (x[[1]]$prior$priorK == "fixed"){
    cat(paste("     ", x[[1]]$prior$Kmax, " component normal mixture estimated using MCMC\n", sep=""))
    cat("     ================================================\n")
  }else{
    cat(paste("     Normal mixture with at most ", x[[1]]$prior$Kmax, " components estimated using RJ-MCMC\n", sep=""))
    cat("     =================================================================")
  }  

  if (ped){
    cat("\nPenalized expected deviance:")
    cat("\n----------------------------\n")
    print(x$PED, ...)
  }  

  if (x[[1]]$prior$priorK != "fixed"){
    cat("\nPosterior distribution of K:")
    cat("\n----------------------------\n")
    Kavail <- unique(c(names(x[[1]]$propK), names(x[[2]]$propK)))
    xpropK <- matrix(0, nrow=2, ncol=length(Kavail))
    rownames(xpropK) <- chNames
    colnames(xpropK) <- Kavail
    xpropK[1, names(x[[1]]$propK)] <- x[[1]]$propK
    xpropK[2, names(x[[2]]$propK)] <- x[[2]]$propK        
    print(xpropK, ...)
  }

  if (dic){
    cat("\nDeviance information criteria:")
    cat("\n------------------------------\n")
    xDIC <- rbind(x[[1]]$DIC, x[[2]]$DIC)
    rownames(xDIC) <- chNames
    print(xDIC, ...)
  }  

  cat("\nPosterior summary statistics for moments of mixture for original data:")
  cat("\n----------------------------------------------------------------------")
  if (x[[1]]$dim == 1){
    if (x[[1]]$nx_w == 1){        
      cat("\nMean:\n")
      xsummyMean <- rbind(x[[1]]$summ.y.Mean, x[[2]]$summ.y.Mean)
      rownames(xsummyMean) <- chNames    
      print(xsummyMean, ...)

      cat("\nStandard deviation:\n")
      xsummySDCorr <- rbind(x[[1]]$summ.y.SDCorr, x[[2]]$summ.y.SDCorr)
      rownames(xsummySDCorr) <- chNames
      print(xsummySDCorr, ...)
    }else{
      cat("\nMeans (chain 1):\n")
      print(x[[1]]$summ.y.Mean, ...)
      cat("\nMeans (chain 2):\n")
      print(x[[2]]$summ.y.Mean, ...)
    
      cat("\nStandard deviations (chain 1):\n")
      print(x[[1]]$summ.y.SDCorr, ...)
      cat("\nStandard deviations (chain 2):\n")
      print(x[[2]]$summ.y.SDCorr, ...)    
    }   
  }else{
    cat("\nMeans (chain 1):\n")
    print(x[[1]]$summ.y.Mean, ...)
    cat("\nMeans (chain 2):\n")
    print(x[[2]]$summ.y.Mean, ...)
    
    cat("\nStandard deviations and correlations (chain 1):\n")
    print(x[[1]]$summ.y.SDCorr, ...)
    cat("\nStandard deviations and correlations (chain 2):\n")
    print(x[[2]]$summ.y.SDCorr, ...)    
  }

  if (!is.null(x[[1]]$summ.expy.Mean)){
    cat("\nPosterior summary statistics for mean of exp(data):")
    cat("\n---------------------------------------------------\n")
    if (x[[1]]$dim == 1 & x[[1]]$nx_w == 1){
      xsummexpyMean <- rbind(x[[1]]$summ.expy.Mean, x[[2]]$summ.expy.Mean)
      rownames(xsummexpyMean) <- chNames    
      print(xsummexpyMean, ...)
    }else{
      cat("\nChain 1:\n")
      print(x[[1]]$summ.expy.Mean, ...)
      cat("\nChain 2:\n")
      print(x[[2]]$summ.expy.Mean, ...)
    }  
  }
    
  return(invisible(x))  
}
