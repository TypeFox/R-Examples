#' Utility functions for working with Random Number Generators
#' 
#' This package contains a set of functions for working with
#' Random Number Generators (RNGs). In particular, it defines a generic
#' S4 framework for getting/setting the current RNG, or RNG data
#' that are embedded into objects for reproducibility.
#' 
#' Notably, convenient default methods greatly facilitate the way current
#' RNG settings can be changed.
#' 
#' @name rngtools
#' @docType package
#' 
#' @import stringr
#' @import digest
#' @import methods
#' @import pkgmaker
#' 
#' @examples 
#' 
#' showRNG()
#' s <- getRNG()
#' RNGstr(s)
#' RNGtype(s)
#' 
#' # get what would be the RNG seed after set.seed
#' s <- nextRNG(1234)
#' showRNG(s)
#' showRNG( nextRNG(1234, ndraw=10) )
#' 
#' # change of RNG kind
#' showRNG()
#' k <- RNGkind()
#' k[2L] <- 'Ahrens'
#' try( RNGkind(k) )
#' setRNG(k)
#' showRNG()
#' # set encoded kind 
#' setRNG(501L)
#' showRNG()
#' 
#' # use as set seed
#' setRNG(1234)
#' showRNG()
#' r <- getRNG()
#' 
#' # extract embedded RNG specifications
#' runif(10)
#' setRNG(list(1, rng=1234))
#' rng.equal(r)
#' 
#' # restore default RNG (e.g., after errors)
#' RNGrecovery()
#' showRNG()
#' 
NULL
