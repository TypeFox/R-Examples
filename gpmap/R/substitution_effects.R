## This file contains functions using allele substitution effects for measuring 
## monotonicity of gp maps.

degree_of_monotonicity <- function(gpmap) {
  # Calculate the degree of monotonicity of gpmap object

  if (!inherits(gpmap,"gpmap")) {
    stop("Input gpmap is not an object of class \"gpmap\". Contruct by calling generate_gpmap.\n")
  }
  nmaps <- gpmap$nmaps
  monodegree <- NULL
  if (nmaps == 1) {
    monodegree <- degree_of_monotonicity_single(gpmap)
  }
  else {
    multidm <- foreach (i = 1:nmaps) %do% { 
      degree_of_monotonicity_single(generate_gpmap(gpmap$values[, i]))
    }
    ## TODO rewrite to use plyr
    monodegree$m.k <- foreach (i = 1:nmaps, .combine = "cbind") %do% { multidm[[i]]$m.k }
    monodegree$s.k.w <- foreach (i = 1:nmaps, .combine = "cbind") %do% { multidm[[i]]$s.k.w }
    monodegree$dm <- foreach (i = 1:nmaps, .combine = "cbind") %do% { multidm[[i]]$dm }
  }
  colnames(monodegree$m.k)      <- gpmap$mapnames
  rownames(monodegree$m.k)      <- gpmap$locinames
  colnames(monodegree$s.k.w)    <- gpmap$mapnames
  rownames(monodegree$s.k.w)    <- gpmap$locinames
  colnames(monodegree$dm)       <- gpmap$mapnames

  gpmap$degree.monotonicity        <- monodegree$dm
  gpmap$degree.monotonicity.locus  <- monodegree$m.k 
  gpmap$locus.weigth               <- monodegree$s.k.w 
  return(gpmap)
}

##Aimed at internal use.
degree_of_monotonicity_single <- function(gpmap) {

  nloci <- gpmap$nloci
  values <- gpmap$values
  mk <- rep(NA, nloci) 
  muk <- rep(NA, nloci)  #new measure linear in Pk and Nk
  s.k.total <- mk
  s.k.w <- mk
  dm <- NA
  
  if (!is.nan(sum(values))) { 
  # reshape GPmap to 3x3x...x3 array with genotype indexes for nloci 
  # biallelic loci
    dim(values) <- rep(3, nloci)
    for (locus in 1:nloci) {
      if (nloci == 1) { 
        subst.effect = diff(values) 
      }
      else {
        g <- 1:3 ^ nloci 
        dim(g) <- rep(3, nloci) 
        sortg <- values[aaply(g, locus, sort)]
        dim(sortg) <- c(3, 3 ^ (nloci - 1))  
        subst.effect <- diff(sortg) 
      }
      # sum of positive substitution effect across backgrounds
      s.plus <- sum(subst.effect[subst.effect>0])   
      # sum of negative substitution effect across backgrounds   	    
      s.minus <- abs(sum(subst.effect[subst.effect<0]))	    
      s.min <- min(c(s.plus, s.minus))
      s.max <- max(c(s.plus, s.minus))
      s.k.total[locus] <- s.plus + s.minus
      if (s.max == 0) { 
        #mk[locus] <- -1 
	mk[locus] <- -1
      }
      else { 
        #mk[locus] <- (s.max - s.min) / s.max 
	mk[locus] <- (s.max - s.min) / (s.max + s.min)
      }
    }
    if (max(mk) == -1) { 
      s.k.total <- rep(1 / nloci, nloci) 
    }
  }
  s.k.w <- s.k.total / sum(s.k.total)
  return(list(#"m.k"        = matrix(mk, ncol = 1),
	      "m.k"        = matrix(mk, ncol = 1),
              "s.k.weight" = matrix(s.k.w, ncol = 1),
	      #"dm"         = s.k.w %*% mk,
	      "dm"         = s.k.w %*% mk))
}
