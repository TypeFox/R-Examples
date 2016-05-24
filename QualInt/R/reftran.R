reftran <- function(refmean = NULL, refstd = NULL, 
                    refrisk = NULL, refrate = NULL,
                    effect, scale = c("RD", "RR", "OR"),
                    type = c("continuous", "binary", "survival")){
  
  type <- match.arg(type)
  
  effect <- as.vector(effect)
  if(!is.numeric(effect))
    stop("effect should be a numeric vector")
  neff <- length(effect)
  
  ref <- switch(type, 
                continuous = as.vector(refmean), 
                binary = as.vector(refrisk), 
                survival = as.vector(refrate))
  if(!is.numeric(ref))
    stop("refmean, refrisk, or refrate should be a numeric vector")
  nsbp <- length(ref)
  
  if(neff < nsbp) {
    effect <- c(effect, rep(0, nsbp - neff))
    warning("the length of effect is shorter than the number of subgroups,
              therefore assumed as 0 by default")
  } else if(neff > nsbp) {
    effect <- effect[1 : neff]
    warning("the length of effect is longer than the number of subgroups,
              therefore only the first nsbp elements will be used")
  }
  
  if(type == "continuous") {
    
    refstd <- as.vector(refstd)
    if(!is.numeric(refstd) | (length(refstd) != nsbp))
      stop("refstd should be a numeric vector with the same length as refmean")
    commean <- ref + effect
    comstd <- refstd
    list(mean = list(ref, commean), std = list(refstd, comstd), risk = NULL,
         rate = NULL)
    
  } else if(type == "binary") {
    
    scale = match.arg(scale)
    comrisk <- switch(scale, 
                      RD = ref + effect, RR = ref * effect, 
                      OR = effect * ref / (1 + (effect - 1) * ref))
    list(mean = NULL, std = NULL, risk = list(ref, comrisk), rate = NULL)
     
  } else if(type == "survival") {
    
    comrate <- 1 - (1 - refrate) ^ effect
    list(mean = NULL, std = NULL, risk = NULL, rate = list(ref, comrate))
     
  }
  
}