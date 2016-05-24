##
##  PURPOSE:   Function to provide chains for mixture components
##             being shifted and scaled into the original scale and also
##             being relabeled if this is requested.
##             Also provide chains for derived quantities like standard deviations
##             and correlations.
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   15/07/2013
##
##  FUNCTION:  NMixChainComp.GLMM_MCMC
##             
## ======================================================================

## *************************************************************
## NMixChainComp.GLMM_MCMC
## *************************************************************
NMixChainComp.GLMM_MCMC <- function(x, relabel = TRUE, param = c("w_b", "mu_b", "var_b", "sd_b", "cor_b", "Sigma_b", "Q_b", "Li_b"))
{
  if (x$prior.b$priorK != "fixed") stop("Not implemented when x$prior.b$priorK is not fixed.")
  param <- match.arg(param)
  if (relabel) order <- x$order_b else order <- matrix(rep(1:x$K_b[1], nrow(x$w_b)), ncol=x$K_b[1], byrow=TRUE)

  ##### Mixture weights
  ##### -----------------------------
  if (param == "w_b"){
    if (relabel){
      w <- matrix(nrow = nrow(x$w_b), ncol = x$K_b[1])
      for (k in 1:x$K_b[1]) w[,k] <- x$w_b[cbind(1:nrow(x$w_b), x$order_b[,k])]
      colnames(w) <- colnames(x$w_b)
      return(w)
    }else{
      return(x$w_b)
    }
  }else{

    ##### Mixture means
    ##### ---------------------------
    if (param == "mu_b"){
      if (relabel){
        order <- x$order_b
      }else{
        order <- matrix(rep(1:x$K_b[1], nrow(x$mu_b)), ncol = x$K_b[1], byrow = TRUE)    
      }
      mu <- matrix(nrow = nrow(x$mu_b), ncol = x$dimb * x$K_b[1])
      for (k in 1:x$K_b[1]){
        for (j in 1:x$dimb){
          i <- (k-1)*x$dimb + j
          mu[,i] <- x$scale.b$shift[j] + x$scale.b$scale[j] * x$mu_b[cbind(1:nrow(x$mu_b), (order[,k] - 1)*x$dimb + j)]
        }
      }
      colnames(mu) <- paste("mu", rep(1:x$K_b[1], each = x$dimb), ".", rep(1:x$dimb, x$K_b[1]), sep = "")
      return(mu)
    }else{

      ##### Quantities derived from the mixture covariance matrices
      ##### -----------------------------------------------------------------------
      
      ### Create replicated scale vector needed to multiply columns of Sigma_b (and after inversion of Q_b and Li_b)
      ### from "left" and from "right" (as if full matrices were multiplies)
      ### to get scaled covariance matrices of random effects (or their inversions or Cholesky factors of inversions).
      names(x$scale.b$scale) <- paste("s", 1:x$dimb, sep = "")
      scaleRepL <- x$scale.b$scale
      if (x$dimb > 1){
        for (j in 1:(x$dimb - 1)) scaleRepL <- c(scaleRepL, x$scale.b$scale[-(1:j)])
      }
      scaleRepR <- rep(x$scale.b$scale, x$dimb:1)

      ### Replicate scaleRepL and scaleRepR K-times
      scaleRepL <- rep(scaleRepL, x$K_b[1])
      scaleRepR <- rep(scaleRepR, x$K_b[1])      

      ### Get rescaled sample
      if (param %in% c("var_b", "sd_b", "cor_b", "Sigma_b")){
        Samp <- matrix(rep(scaleRepL, nrow(x$Sigma_b)), nrow = nrow(x$Sigma_b), byrow = TRUE) * x$Sigma_b * matrix(rep(scaleRepR, nrow(x$Sigma_b)), nrow = nrow(x$Sigma_b), byrow = TRUE)
      }else{
        switch (param,
          "Q_b"  = {Samp <- matrix(rep(1 / scaleRepL, nrow(x$Q_b)), nrow = nrow(x$Q_b), byrow = TRUE) * x$Q_b * matrix(rep(1 / scaleRepR, nrow(x$Q_b)), nrow = nrow(x$Q_b), byrow = TRUE)},
          "Li_b" = {Samp <- matrix(rep(1 / scaleRepL, nrow(x$Li_b)), nrow = nrow(x$Li_b), byrow = TRUE) * x$Li_b}
        )
      }

      ### Get relabeled sample (if needed)
      LTp <- (x$dimb * (x$dimb + 1)) / 2
      
      if (relabel){        
        SampRel <- matrix(nrow = nrow(Samp), ncol = ncol(Samp))
        colnames(SampRel) <- colnames(Samp)
        for (k in 1:x$K_b[1]){
          for (j in 1:LTp){
            i <- (k - 1) * LTp + j
            SampRel[, i] <- Samp[cbind(1:nrow(Samp), (x$order_b[,k] - 1) * LTp + j)]
          }  
        }
        Samp <- SampRel
        rm(list = "SampRel")        
      }

      ### Return what is requested
      Idiag <- matrix(0, nrow = x$dimb, ncol = x$dimb)
      Idiag[lower.tri(Idiag, diag = TRUE)] <- 1:LTp
      jdiag <- diag(Idiag)                                                  ## indeces of diagonal elements in a lower triangle
      jdiagAll <- numeric()
      for (k in 1:x$K_b[1]) jdiagAll <- c(jdiagAll, (k - 1)*LTp + jdiag)    ## column indeces in Samp corresponding to diagonal elements
      
      if (param %in% c("Sigma_b", "Q_b", "Li_b")){
        return(Samp)
      }else{
        switch (param,
          "var_b" = {
            Samp <- Samp[, jdiagAll]
            colnames(Samp) <- paste("var", rep(1:x$K_b[1], each = x$dimb), ".", rep(1:x$dimb, x$K_b[1]), sep = "")
            return(Samp)
          },
          "sd_b"  = {
            Samp <- sqrt(Samp[, jdiagAll])
            colnames(Samp) <- paste("sd", rep(1:x$K_b[1], each = x$dimb), ".", rep(1:x$dimb, x$K_b[1]), sep = "")
            return(Samp)
          },
          "cor_b" = {
            if (x$dimb == 1) return(numeric(0))
            
            ncor <- (x$dimb - 1) * x$dimb / 2
            cSamp <- matrix(nrow = nrow(Samp), ncol = ncor * x$K_b[1])
            cName <- character(ncor * x$K_b[1])
            cc <- 0
            for (k in 1:x$K_b[1]){
              for (j in 1:(x$dimb - 1)){
                for (i in (j+1):x$dimb){
                  cc <- cc + 1
                  cName[cc] <- paste("cor", k, ".", i, ".", j, sep = "")
                  cSamp[, cc] <- Samp[, paste("Sigma", k, ".", i, ".", j, sep = "")] / sqrt(Samp[, paste("Sigma", k, ".", i, ".", i, sep = "")] * Samp[, paste("Sigma", k, ".", j, ".", j, sep = "")])
                }  
              }  
            }
            colnames(cSamp) <- cName
            return(cSamp)
          }
        )
      }        
    }  
  }    
}
