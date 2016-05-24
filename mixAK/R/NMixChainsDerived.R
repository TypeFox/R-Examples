##
##  PURPOSE:   Add derived chains to the object computed using the NMixMCMC function
##             + basic posterior summary statistics for them
##
##  AUTHOR:    Arnost Komarek (LaTeX: Arno\v{s}t Kom\'arek)
##             arnost.komarek[AT]mff.cuni.cz
##
##  CREATED:   09/06/2009
##
##  FUNCTIONS:  NMixChainsDerived
##
## ======================================================================

## *************************************************************
## NMixChainsDerived
## *************************************************************
NMixChainsDerived <- function(object)
{
  thispackage <- "mixAK"

  if (!(class(object) %in% c("NMixMCMClist", "NMixMCMC"))) stop("object must be of class NMixMCMClist or NMixMCMC")
  
  ##### Quantiles and names of posterior summary statistics which are computed
  ##### for derived quantities
  qProbs <- c(0, 0.025, 0.25, 0.5, 0.75, 0.975, 1)
  nSumm <- c("Mean", "Std.Dev.", "Min.", "2.5%", "1st Qu.", "Median", "3rd Qu.", "97.5%", "Max.")

  ##### Function which works out one chain
  NMixChDer <- function(obj1)
  {
    LTp <- obj1$dim * (obj1$dim + 1)/2
    Krandom <- 1*(obj1$prior$priorK != "fixed")
    if (Krandom){
      M <- length(obj1$K)
      RES <- matrix(
               .C("NMix_ChainsDerived",
                       chEexpY    = double(M * obj1$dim),
                       dwork      = double(LTp),
                       err        = integer(1),
                       p          = as.integer(obj1$dim),
                       shiftScale = as.double(c(obj1$scale$shift, obj1$scale$scale)),
                       chK        = as.integer(obj1$K),
                       chw        = as.double(obj1$w),
                       chmu       = as.double(obj1$mu),
                       chLi       = as.double(obj1$Li),
                       M          = as.integer(M),
                       Krandom    = as.integer(Krandom),
                 PACKAGE=thispackage)$chEexpY,
                 ncol = obj1$dim, byrow=TRUE)
    }else{      
      M <- length(obj1$w) / obj1$K[1]
      
      RES <- matrix(
               .C("NMix_ChainsDerived",
                      chEexpY    = double(M * obj1$dim),
                      dwork      = double(LTp),
                      err        = integer(1),
                      p          = as.integer(obj1$dim),
                      shiftScale = as.double(c(obj1$scale$shift, obj1$scale$scale)),
                      chK        = as.integer(obj1$K[1]),
                      chw        = as.double(t(obj1$w)),
                      chmu       = as.double(t(obj1$mu)),
                      chLi       = as.double(t(obj1$Li)),
                      M          = as.integer(M),
                      Krandom    = as.integer(Krandom),
                PACKAGE=thispackage)$chEexpY,
                ncol = obj1$dim, byrow=TRUE)
    }
    colnames(RES) <- paste("expy.Mean.", 1:obj1$dim, sep="")
    obj1$chains.derived <- as.data.frame(RES)

    if (obj1$dim == 1){
      PMean  <- mean(obj1$chains.derived[, "expy.Mean.1"], na.rm=TRUE)
      PSD    <- sd(obj1$chains.derived[, "expy.Mean.1"], na.rm=TRUE)
      PQuant <- quantile(obj1$chains.derived[, "expy.Mean.1"], prob=qProbs, na.rm=TRUE)
      obj1$summ.expy.Mean <- c(PMean, PSD, PQuant)
      names(obj1$summ.expy.Mean) <- nSumm
    }else{
      Naam <- paste("expy.Mean.", 1:obj1$dim, sep="")
      PMean <- apply(obj1$chains.derived[, Naam], 2, mean, na.rm=TRUE)
      PSD <- apply(obj1$chains.derived[, Naam], 2, sd, na.rm=TRUE)
      PQuant <- apply(obj1$chains.derived[, Naam], 2, quantile, prob=qProbs, na.rm=TRUE)      
      obj1$summ.expy.Mean <- rbind(PMean, PSD, PQuant)
      rownames(obj1$summ.expy.Mean) <- nSumm      
    }

    return(obj1)
  }  

  
  ##### NMixMCMClist
  if (class(object) == "NMixMCMClist"){
    if (object[[1]]$nx_w > 1) stop("This function has not (yet) been implemented if a factor covariate on mixture weights is present.")
    for (ch in 1:2){
      object[[ch]] <- NMixChDer(object[[ch]])      
    }
    return(object)
  }  


  ##### NMixMCMC
  if (class(object) == "NMixMCMC"){
    if (object$nx_w > 1) stop("This function has not (yet) been implemented if a factor covariate on mixture weights is present.")      
    return(NMixChDer(object))
  }
}  
