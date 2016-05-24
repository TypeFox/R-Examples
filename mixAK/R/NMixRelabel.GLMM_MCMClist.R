##
##  PURPOSE:   Re-labeling of the MCMC output.
##             * method for objects of class GLMM_MCMClist
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   25/06/2013
##             16/08/2013  'jointly' option added
##
##  FUNCTION:  NMixRelabel.GLMM_MCMClist (25/06/2013) 
##
## ======================================================================

## *************************************************************
## NMixRelabel.GLMM_MCMClist
## *************************************************************
NMixRelabel.GLMM_MCMClist <- function(object, type = c("mean", "weight", "stephens"), par,
                                      prob = c(0.025, 0.5, 0.975), keep.comp.prob = FALSE,
                                      jointly = FALSE, info, silent = FALSE,
                                      parallel = FALSE, ...)
{
  if (jointly){
    keepMCMC1 <- length(object[[1]]$w_b) / object[[1]]$K_b[1]
    keepMCMC2 <- length(object[[2]]$w_b) / object[[2]]$K_b[1]    
    
    objBoth <- object[[1]]
    l_alpha <- sum(objBoth$p) + sum(objBoth$fixed.intercept)
    
    if (is.null(object[[1]]$K_b) | is.null(object[[2]]$K_b) | is.null(object[[1]]$w_b) | is.null(object[[2]]$w_b) |
        is.null(object[[1]]$mu_b) | is.null(object[[2]]$mu_b) | is.null(object[[1]]$Li_b) | is.null(object[[2]]$Li_b) |
        is.null(object[[1]]$Q_b) | is.null(object[[2]]$Q_b) | is.null(object[[1]]$Sigma_b) | is.null(object[[1]]$Sigma_b)){
      stop("object does not contain sampled values of mixture related parameter(s)")
    }
    objBoth$K_b     <- c(objBoth$K_b, object[[2]]$K_b)
    objBoth$w_b     <- c(objBoth$w_b, object[[2]]$w_b)
    objBoth$mu_b    <- rbind(objBoth$mu_b,    object[[2]]$mu_b)
    objBoth$Li_b    <- rbind(objBoth$Li_b,    object[[2]]$Li_b)
    objBoth$Q_b     <- rbind(objBoth$Q_b,     object[[2]]$Q_b)
    objBoth$Sigma_b <- rbind(objBoth$Sigma_b, object[[2]]$Sigma_b)    
    
    if (objBoth$R["Rc"]){
      if (is.null(object[[1]]$sigma_eps) | is.null(object[[2]]$sigma_eps)) stop("object does not contain sampled values of sigma(eps)")
      objBoth$sigma_eps <- rbind(objBoth$sigma_eps, object[[2]]$sigma_eps)      
    }
    
    if (l_alpha){
      if (is.null(object[[1]]$alpha) | is.null(object[[2]]$alpha)) stop("object does not contain sampled values of fixed effects")
      objBoth$alpha <- rbind(objBoth$alpha, object[[2]]$alpha)              
    }

    objBoth <- NMixRelabel.GLMM_MCMC(objBoth, type = type, par = par, prob = prob, keep.comp.prob = keep.comp.prob, info = info, silent = silent, ...)

    object[[1]]$order_b <- objBoth$order_b[1:keepMCMC1,]
    object[[2]]$order_b <- objBoth$order_b[(keepMCMC1 + 1):(keepMCMC1 + keepMCMC2),]

    object[[1]]$rank_b <- objBoth$rank_b[1:keepMCMC1,]
    object[[2]]$rank_b <- objBoth$rank_b[(keepMCMC1 + 1):(keepMCMC1 + keepMCMC2),]

    object[[1]]$poster.comp.prob_u <- object[[2]]$poster.comp.prob_u <- objBoth$poster.comp.prob_u
    object[[1]]$poster.comp.prob_b <- object[[2]]$poster.comp.prob_b <- objBoth$poster.comp.prob_b
    object[[1]]$poster.comp.prob   <- object[[2]]$poster.comp.prob   <- objBoth$poster.comp.prob

    if (length(prob)){
      object[[1]]$quant.comp.prob_b <- object[[2]]$quant.comp.prob_b <- list()
      object[[1]]$quant.comp.prob   <- object[[2]]$quant.comp.prob   <- list()      
      for (i in 1:length(prob)){
        object[[1]]$quant.comp.prob_b[[i]] <- object[[2]]$quant.comp.prob_b[[i]] <- objBoth$quant.comp.prob_b[[i]]
        object[[1]]$quant.comp.prob[[i]]   <- object[[2]]$quant.comp.prob[[i]]   <- objBoth$quant.comp.prob[[i]]
      }
      names(object[[1]]$quant.comp.prob_b) <- names(object[[2]]$quant.comp.prob_b) <- names(objBoth$quant.comp.prob_b)
      names(object[[1]]$quant.comp.prob)   <- names(object[[2]]$quant.comp.prob)   <- names(objBoth$quant.comp.prob)
     }else{
      object[[1]]$quant.comp.prob_b <- object[[2]]$quant.comp.prob_b <- NULL
      object[[1]]$quant.comp.prob   <- object[[2]]$quant.comp.prob   <- NULL      
    }  
    
    if (keep.comp.prob){
      object[[1]]$comp.prob_b <- objBoth$comp.prob_b[1:keepMCMC1,]
      object[[2]]$comp.prob_b <- objBoth$comp.prob_b[(keepMCMC1 + 1):(keepMCMC1 + keepMCMC2),]

      object[[1]]$comp.prob <- objBoth$comp.prob[1:keepMCMC1,]
      object[[2]]$comp.prob <- objBoth$comp.prob[(keepMCMC1 + 1):(keepMCMC1 + keepMCMC2),]            
    }else{
      object[[1]]$comp.prob_b <- object[[2]]$comp.prob_b <- NULL
      object[[1]]$comp.prob   <- object[[2]]$comp.prob   <- NULL      
    }          

    object[[1]]$poster.mean.w_b     <- object[[2]]$poster.mean.w_b     <- objBoth$poster.mean.w_b
    object[[1]]$poster.mean.mu_b    <- object[[2]]$poster.mean.mu_b    <- objBoth$poster.mean.mu_b
    object[[1]]$poster.mean.Q_b     <- object[[2]]$poster.mean.Q_b     <- objBoth$poster.mean.Q_b
    object[[1]]$poster.mean.Sigma_b <- object[[2]]$poster.mean.Sigma_b <- objBoth$poster.mean.Sigma_b
    object[[1]]$poster.mean.Li_b    <- object[[2]]$poster.mean.Li_b    <- objBoth$poster.mean.Li_b    
        
    return(object)
  }

  ##### The remaining code applies only if jointly = FALSE
  if (parallel){
    RAlg <- NMixRelabelAlgorithm(type = type, par = par, dim = object[[1]]$dimb)
    if (missing(info)) info <- object[[1]]$nMCMC["keep"]
    
    #require("parallel")
    if (parallel::detectCores() < 2) warning("It does not seem that at least 2 CPU cores are available needed for efficient parallel re-labelling of the two chains.")      
    cl <- parallel::makeCluster(2)      
    if (!silent){
      cat("\nParallel re-labelling of the two chains\n")
      cat("=======================================\n\n")
    }      
    tmpObj <- list(object[[1]], object[[2]])
    class(tmpObj) <- "GLMM_MCMClist"
    tmpObj <- parallel::parLapply(cl, tmpObj, NMixRelabel.GLMM_MCMC,
                                  type = RAlg$relabel$type, par = RAlg$relabel$par, prob = prob, keep.comp.prob = keep.comp.prob, info = info, silent = silent, ...)        
    parallel::stopCluster(cl)

    elemObj <- names(object)[-(1:2)]
    for (i in 1:length(elemObj)) tmpObj[[elemObj[i]]] <- object[[elemObj[i]]]

    class(tmpObj) <- "GLMM_MCMClist"    
    return(tmpObj)
  }else{
    if (!silent){
      cat("\nRe-labelling chain number 1\n")
      cat("===========================\n")
    }  
    object[[1]] <- NMixRelabel(object[[1]], type = type, par = par, prob = prob, keep.comp.prob = keep.comp.prob, info = info, silent = silent, ...)

    if (!silent){    
      cat("\nRe-labelling chain number 2\n")
      cat("===========================\n")
    }  
    object[[2]] <- NMixRelabel(object[[2]], type = type, par = par, prob = prob, keep.comp.prob = keep.comp.prob, info = info, silent = silent, ...)
    cat("\n")
    
    return(object)
  }
}  
