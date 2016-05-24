## This file contains functions using isotone regression for measuring 
## monotonicity in gpmap objects.

decompose_monotone <- function(gpmap) {
    
  if (!inherits(gpmap, "gpmap")) { 
    stop("Create gpmap by calling the function generategpmap first. \n") 
  }
  nloci <- gpmap$nloci
  nmaps <- gpmap$nmaps
  dmap <- NULL
  if (nmaps == 1) {
    dmap <- NULL
    if (nloci > 3) {
      dmap <- decompose_monotone_single(gpmap)
    }
    else {
      dmap <- decompose_monotone_single(gpmap)
    }
  }
  else {
    tmp <- foreach (i = 1:nmaps) %dopar% {
      decompose_monotone_single(generate_gpmap(gpmap$values[, i]))
    }
    dmap$monoR2 <- foreach (m = 1:nmaps, .combine = "cbind") %do% {
      matrix(tmp[[m]]$monoR2)
    }
    dmap$values.mono <- foreach (d = 1:nmaps, .combine = "cbind") %do% {
      tmp[[d]]$values.mono
    }
  }
  colnames(dmap$values.mono) <- gpmap$mapnames
  colnames(dmap$monoR2)   <- gpmap$mapnames

  gpmap$monoR2 <- dmap$monoR2          
  gpmap$values.mono <- dmap$values.mono
  return(gpmap)
}


## Low-level functions ##
decompose_monotone_single <- function(gpmap) {
  
  nloci <- gpmap$nloci

  ## construct all possible combinations of plusalleles
  tmpargs <- NULL 
  for (i in 1:nloci) {
    tmpargs[[paste('Locus', i, sep = "")]] <- c(1, 2)
  }
  plusallele <- do.call(expand.grid, tmpargs)
  decomposed.all  <- foreach (i = 1:2^nloci) %dopar% {
    monotone_regression(gpmap, plusallele[i, ])
  }   
  decomposed <- decomposed.all[[which.min(laply(decomposed.all, 
                                                function(x) return(x$fval)))]]
  return(list("values.mono" = matrix(decomposed$x, ncol = 1),
              "monoR2"  = var(decomposed$x) / var(decomposed$y)))
}

### Low-level functions ### 

monotone_regression <- function(gpmap, plusallele) {
  # Wrapper of the monotone regression functions. 

  ## build up partial ordering of genotype space base on input vector 
  ## plusallele
  nloci <- length(plusallele)
  apart <- partial_genotype_order(plusallele)

  ## do monotone_regression of genotypic values using isotone::activeSet
  fit.ls <- activeSet(apart, "LS", y = gpmap$values, weights = rep(1, 3 ^ nloci))
  return(fit.ls)
}

#treat ´m´,´d´ and ´i´ as global during checking (to avoid NOTE)
if(getRversion() >= "2.15.1")  utils::globalVariables(c("i","m","d"))
