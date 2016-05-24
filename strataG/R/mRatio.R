#' @title M ratio
#' @description Calculate Garza-Williamson M ratio (bottleneck) statistic for 
#'   microsattelite data.
#'
#' @param g a \linkS4class{gtypes} object.
#' @param by.strata calculate ratio for each stratum separately?
#' @param rpt.size set of values to check for allele repeat size. Function 
#'   will use the largest common denominator found in this vector or return \code{NA}.
#'
#' @return If \code{by.strata = TRUE}, a matrix with loci as rows and strata 
#'   as columns, otherwise a vector with the statistic for each locus over all
#'   samples.
#'
#' @note The function will only compute the metric for microastellite loci, which 
#'  is defined as loci with allele labels that can be converted to numeric 
#'  values in their entirety and have a fixed repeat size. \code{NA} is returned 
#'  for all others.
#'
#' @references Garza, J.C. and E.G. Williamson. 2001. Detection of reduction in 
#'   population size using data form microsatellite loci. 
#'   Molecular Ecology 10(2):305-318.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} (adapted from code by Sean Hoban)
#'
#' @examples
#' data(msats.g)
#' 
#' m.by.strata <- mRatio(msats.g, TRUE)
#' m.by.strata
#' 
#' m.overall <- mRatio(msats.g, FALSE)
#' m.overall
#' 
#' @importFrom Hmisc all.is.numeric
#' @export

mRatio <- function(g, by.strata = TRUE, rpt.size = 8:2) {
  # function takes a matrix of allele frequencies (can be more than one column)
  # rownames represent allele sizes
  calc.mratio<- function(freqs) {
    # extract first column
    freqs <- freqs[, 1]
    
    if(length(freqs) == 1) {
      warning("only one allele")
      return(NA)
    }
    
    # check if rownames are numerics
    if(!all.is.numeric(names(freqs))) {
      warning("allele names are non-numeric")
      return(NA)
    }
    
    if(all(freqs == 0)) { 
      warning("all frequencies are 0")
      NA
    } else {
      # sort alleles in numerical order
      freqs <- freqs[order(as.numeric(names(freqs)))]
      # convert names to numbers to get sizes
      sizes <- as.numeric(names(freqs))
      # find repeat sizes
      size.diff <- diff(sizes)
      rpt.found <- FALSE
      for(r in sort(rpt.size, decreasing = TRUE)) {
        if(all(size.diff %% r == 0)) {
          rpt.found <- TRUE
          break
        }
      }
      if(!rpt.found) {
        warning("valid repeat length not found")
        return(NA)
      }
      # find smallest and largest alleles that are present
      smallest <- sizes[which.min(sizes[freqs > 0])]
      largest <- sizes[which.max(sizes[freqs > 0])]
      # compute number of alleles between smallest and largest
      n <- (largest - smallest) / r
      # calculate metric
      sum(freqs > 0) / (n + 1)
    }
  }
  
  if(nStrata(g) == 1 & by.strata) {
    by.strata <- FALSE
    g <- g[, , strataNames(g)]
  }
  if(by.strata) {
    freqs <- alleleFreqs(g, by.strata = TRUE)
    do.call(rbind, lapply(freqs, function(loc) apply(loc, 3, calc.mratio)))
  } else {
    freqs <- alleleFreqs(g, by.strata = FALSE)
    sapply(freqs, calc.mratio)
  }
}