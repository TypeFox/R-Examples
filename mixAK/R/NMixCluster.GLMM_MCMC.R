##
##  PURPOSE:   Clustering based on the MCMC output.
##             * method for objects of class GLMM_MCMC
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   20/08/2014
##
##  FUNCTION:  NMixCluster.GLMM_MCMC (20/08/2014) 
##
## ======================================================================

## *************************************************************
## NMixCluster.GLMM_MCMC
## *************************************************************
NMixCluster.GLMM_MCMC <- function(object,
                                  prob = c("poster.comp.prob", "quant.comp.prob", "poster.comp.prob_b", "quant.comp.prob_b", "poster.comp.prob_u"),
                                  pquant = 0.5,
                                  HPD = FALSE,
                                  pHPD = 0.95,
                                  pthresh = -1,
                                  unclass.na = FALSE, ...)
{
  if (object[["prior.b"]][["priorK"]] != "fixed") stop("not implemented when the number of mixture components was not fixed.")
  K <- object[["prior.b"]][["Kmax"]]
  nG <- ifelse(pthresh > 0, ifelse(unclass.na, K, K + 1), K)
  
  prob <- match.arg(prob)
  if (prob == "poster.comp.prob_u") HPD <- FALSE
  if (prob %in% c("poster.comp.prob", "quant.comp.prob") & is.null(object[["comp.prob"]])) HPD <- FALSE
  if (prob %in% c("poster.comp.prob_b", "quant.comp.prob_b") & is.null(object[["comp.prob_b"]])) HPD <- FALSE  
    
  if (substr(prob, 1, 5) == "quant"){
    phat <- object[[prob]][[paste(pquant*100, "%", sep = "")]]
  }else{
    phat <- object[[prob]]
  }    

  if (is.null(phat)) stop("component probabilities not available within the supplied object.")

  pGroup <- apply(phat, 1, max)
  Group <- apply(phat, 1, which.max)
 
  if (HPD){
    if (prob %in% c("poster.comp.prob", "quant.comp.prob")){
      hpd <- coda::HPDinterval(coda::mcmc(object[["comp.prob"]]), prob = pHPD)
    }
    if (prob %in% c("poster.comp.prob_b", "quant.comp.prob_b")){
      hpd <- coda::HPDinterval(coda::mcmc(object[["comp.prob_b"]]), prob = pHPD)
    }

    hpd.low <- matrix(hpd[, "lower"], ncol = K, byrow = TRUE)
    hpd.upp <- matrix(hpd[, "upper"], ncol = K, byrow = TRUE)    
    ihpd.low <- hpd.low[cbind(1:nrow(hpd.low), Group)]
    ihpd.upp <- hpd.upp[cbind(1:nrow(hpd.upp), Group)]    

    Group[ihpd.low < pthresh] <- ifelse(unclass.na, NA, K + 1)    
    RET <- data.frame(Group = Group, fGroup = factor(Group, levels = 1:nG), p = pGroup, p.lower = ihpd.low, p.upper = ihpd.upp)    
  }else{
    Group[pGroup < pthresh] <- ifelse(unclass.na, NA, K + 1)
    RET <- data.frame(Group = Group, fGroup = factor(Group, levels = 1:nG), p = pGroup)
  }    

  return(RET)
}
