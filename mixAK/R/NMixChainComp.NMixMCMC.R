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
##  CREATED:   16/07/2013
##             27/03/2015  mild revision to allow for factor covariates on mixture weights
##
##  FUNCTION:  NMixChainComp.NMixMCMC
##             
## ======================================================================

## *************************************************************
## NMixChainComp.NMixMCMC
## *************************************************************
NMixChainComp.NMixMCMC <- function(x, relabel = TRUE, param = c("w", "mu", "var", "sd", "cor", "Sigma", "Q", "Li"))
{
  if (x$prior$priorK != "fixed") stop("Not implemented when x$prior$priorK is not fixed.")
  param <- match.arg(param)
  if (relabel) order <- x$order else order <- matrix(rep(1:x$K[1], nrow(x$w)), ncol=x$K[1], byrow=TRUE)

  ##### Mixture weights
  ##### -----------------------------
  if (param == "w"){
    if (relabel){
      if (x$nx_w == 1){        
        w <- matrix(nrow = nrow(x$w), ncol = x$K[1])
        for (k in 1:x$K[1]) w[,k] <- x$w[cbind(1:nrow(x$w), x$order[,k])]
        colnames(w) <- colnames(x$w)
        return(w)
      }else{
        w <- matrix(nrow = nrow(x$w), ncol = x$K[1] * x$nx_w)
        for (ixw in 1:x$nx_w){
          wixw <- x$w[, (ixw-1)*x$K[1] + (1:x$K[1])]
          for (k in 1:x$K[1]) w[, (ixw-1)*x$K[1] + k] <- wixw[cbind(1:nrow(wixw), x$order[,k])]          
          rm(list = "wixw")
        }
        colnames(w) <- colnames(x$w)
        return(w)
      }    
    }else{
      return(x$w)
    }
  }else{

    ##### Mixture means
    ##### ---------------------------
    if (param == "mu"){
      if (relabel){
        order <- x$order
      }else{
        order <- matrix(rep(1:x$K[1], nrow(x$mu)), ncol = x$K[1], byrow = TRUE)    
      }
      mu <- matrix(nrow = nrow(x$mu), ncol = x$dim * x$K[1])
      for (k in 1:x$K[1]){
        for (j in 1:x$dim){
          i <- (k-1)*x$dim + j
          mu[,i] <- x$scale$shift[j] + x$scale$scale[j] * x$mu[cbind(1:nrow(x$mu), (order[,k] - 1)*x$dim + j)]
        }
      }
      colnames(mu) <- paste("mu", rep(1:x$K[1], each = x$dim), ".", rep(1:x$dim, x$K[1]), sep = "")
      return(mu)
    }else{

      ##### Quantities derived from the mixture covariance matrices
      ##### -----------------------------------------------------------------------
      
      ### Create replicated scale vector needed to multiply columns of Sigma (and after inversion of Q and Li)
      ### from "left" and from "right" (as if full matrices were multiplies)
      ### to get scaled covariance matrices of random effects (or their inversions or Cholesky factors of inversions).
      names(x$scale$scale) <- paste("s", 1:x$dim, sep = "")
      scaleRepL <- x$scale$scale
      if (x$dim > 1){
        for (j in 1:(x$dim - 1)) scaleRepL <- c(scaleRepL, x$scale$scale[-(1:j)])
      }
      scaleRepR <- rep(x$scale$scale, x$dim:1)

      ### Replicate scaleRepL and scaleRepR K-times
      scaleRepL <- rep(scaleRepL, x$K[1])
      scaleRepR <- rep(scaleRepR, x$K[1])      

      ### Get rescaled sample
      if (param %in% c("var", "sd", "cor", "Sigma")){
        Samp <- matrix(rep(scaleRepL, nrow(x$Sigma)), nrow = nrow(x$Sigma), byrow = TRUE) * x$Sigma * matrix(rep(scaleRepR, nrow(x$Sigma)), nrow = nrow(x$Sigma), byrow = TRUE)
      }else{
        switch (param,
          "Q"  = {Samp <- matrix(rep(1 / scaleRepL, nrow(x$Q)), nrow = nrow(x$Q), byrow = TRUE) * x$Q * matrix(rep(1 / scaleRepR, nrow(x$Q)), nrow = nrow(x$Q), byrow = TRUE)},
          "Li" = {Samp <- matrix(rep(1 / scaleRepL, nrow(x$Li)), nrow = nrow(x$Li), byrow = TRUE) * x$Li}
        )
      }

      ### Get relabeled sample (if needed)
      LTp <- (x$dim * (x$dim + 1)) / 2
      
      if (relabel){        
        SampRel <- matrix(nrow = nrow(Samp), ncol = ncol(Samp))
        colnames(SampRel) <- colnames(Samp)
        for (k in 1:x$K[1]){
          for (j in 1:LTp){
            i <- (k - 1) * LTp + j
            SampRel[, i] <- Samp[cbind(1:nrow(Samp), (x$order[,k] - 1) * LTp + j)]
          }  
        }
        Samp <- SampRel
        rm(list = "SampRel")        
      }

      ### Return what is requested
      Idiag <- matrix(0, nrow = x$dim, ncol = x$dim)
      Idiag[lower.tri(Idiag, diag = TRUE)] <- 1:LTp
      jdiag <- diag(Idiag)                                                  ## indeces of diagonal elements in a lower triangle
      jdiagAll <- numeric()
      for (k in 1:x$K[1]) jdiagAll <- c(jdiagAll, (k - 1)*LTp + jdiag)    ## column indeces in Samp corresponding to diagonal elements
      
      if (param %in% c("Sigma", "Q", "Li")){
        return(Samp)
      }else{
        switch (param,
          "var" = {
            Samp <- Samp[, jdiagAll]
            colnames(Samp) <- paste("var", rep(1:x$K[1], each = x$dim), ".", rep(1:x$dim, x$K[1]), sep = "")
            return(Samp)
          },
          "sd"  = {
            Samp <- sqrt(Samp[, jdiagAll])
            colnames(Samp) <- paste("sd", rep(1:x$K[1], each = x$dim), ".", rep(1:x$dim, x$K[1]), sep = "")
            return(Samp)
          },
          "cor" = {
            if (x$dim == 1) return(numeric(0))
            
            ncor <- (x$dim - 1) * x$dim / 2
            cSamp <- matrix(nrow = nrow(Samp), ncol = ncor * x$K[1])
            cName <- character(ncor * x$K[1])
            cc <- 0
            for (k in 1:x$K[1]){
              for (j in 1:(x$dim - 1)){
                for (i in (j+1):x$dim){
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
