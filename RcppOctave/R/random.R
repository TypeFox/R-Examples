# Function related to Random Number Generators
# 
# Author: Renaud Gaujoux
# Creation: 04 Nov 2011
###############################################################################


#' Drawing from R Uniform Distribution in Octave
#' 
#' This function wraps a call to the standard Octave function \code{rand}, 
#' which is redefined by \code{RcppOctave} to call the R base function 
#' \code{\link[stats]{runif}}.
#' This enables to exactly reproduce stochastic computations in R and Octave, 
#' without changing the original Octave/Matlab code.
#' 
#' @section Difference with plain \code{runif}:
#' Since calling \code{o_runif} or \code{\link[stats]{runif}} is equivalent, this 
#' function may not be really useful for the end user, and is defined for 
#' testing purposes essentially.
#' One possible advantage over plain \code{\link[stats]{runif}} however, is that it 
#' can generate random matrices, instead of only vectors (see examples).
#' 
#' @section Seeding:
#' Because the RNG of R is called used, seeding computations is achieved by a    
#' standard call to \code{\link{set.seed}}.
#' 
#' @section Octave details:
#' \code{RcppOctave} defines a set of functions like \code{rand} that shadow 
#' Octave built-in functions.
#' These functions are defined in the Octave module \code{Rrng.oct} that is 
#' stored in the package \emph{modules/} sub-directory.
#' See \code{Octave.info('modules')} to see this directory's full path.
#' 
#' @param n number of output rows
#' @param p number of output columns (default to \code{n})
#' 
#' @return a numeric vector or a matrix
#' 
#' @templateVar name rand
#' @template OctaveDoc
#' 
#' @family orandom 
#' @export
#' @seealso runif
#' @examples
#' 
#' # Draw random uniform values (in vector form)
#' set.seed(123)
#' o_runif(1)
#' o_runif(1, 10)
#' # The result is identical as calling runif  
#' set.seed(123)
#' runif(1)
#' runif(10)
#' 
#' # Draw random uniform values (in matrix form)
#' set.seed(123)
#' o_runif(2)
#' o_runif(2, 5)
#' 
#' 
o_runif <- function(n, p=n){
	.CallOctave('rand', n, p)	
}

#' Drawing from R Normal Distribution in Octave
#' 
#' This function wraps a call to the standard Octave function \code{randn}, 
#' which is redefined by \code{RcppOctave} to call the R base function 
#' \code{\link[stats]{rnorm}}.
#' This enables to exactly reproduce stochastic computations in R and Octave, 
#' without changing the original Octave/Matlab code.
#' See \code{\link{o_runif}} for more details.
#' 
#' @inheritParams o_runif
#' 
#' @templateVar name randn
#' @template OctaveDoc
#' @export
#' @family orandom
#' @seealso rnorm
#' @examples
#' 
#' # Draw random normal values (in vector form)
#' set.seed(123)
#' o_rnorm(1)
#' o_rnorm(1, 10)
#' 
#' # Draw random normal values (in matrix form)
#' set.seed(123)
#' o_rnorm(2)
#' o_rnorm(2, 5)
#' 
o_rnorm <- function(n, p=n){
	.CallOctave('randn', n, p)
}

#' Drawing from R Exponential Distribution in Octave
#' 
#' This function wraps a call to the standard Octave function \code{rande}, 
#' which is redefined by \code{RcppOctave} to call the R base function 
#' \code{\link[stats]{rexp}}.
#' This enables to exactly reproduce stochastic computations in R and Octave, 
#' without changing the original Octave/Matlab code.
#' See \code{\link{o_runif}} for more details.
#'  
#' @inheritParams o_runif 
#' 
#' @templateVar name rande
#' @template OctaveDoc
#' @export
#' @family orandom
#' @seealso rexp
#' @examples
#' 
#' # Draw random exponential values (in vector form)
#' set.seed(123)
#' o_rexp(1)
#' o_rexp(1, 10)
#' 
#' # Draw random normal values (in matrix form)
#' set.seed(123)
#' o_rexp(2)
#' o_rexp(2, 5)
#' 
o_rexp <- function(n, p=n){
	.CallOctave('rande', n, p)
}

#' Drawing from R Gamma Distribution in Octave
#' 
#' This function wraps a call to the standard Octave function \code{randg}, 
#' which is redefined by \code{RcppOctave} to call the R base function 
#' \code{\link[stats]{rgamma}}.
#' This enables to exactly reproduce stochastic computations in R and Octave, 
#' without changing the original Octave/Matlab code.
#' See \code{\link{o_runif}} for more details.
#' 
#' @inheritParams o_runif 
#' @param shape Mean of the Gamma distribution
#' @param scale Scale of the Gamma distribution
#' 
#' @templateVar name randg
#' @template OctaveDoc
#' 
#' @export
#' @family orandom
#' @seealso rgamma 
#' @examples
#' 
#' # Draw random gamma values (in vector form)
#' set.seed(123)
#' o_rgamma(1)
#' o_rgamma(1, 10)
#' 
#' # Draw random gamma values (in matrix form)
#' set.seed(123)
#' o_rgamma(2)
#' o_rgamma(2, 5)
#' 
#' # Draw random gamma values with shape and scale parameters
#' o_rgamma(1, 5, shape=2)
#' o_rgamma(1, 10, scale=0.5)
#' 
o_rgamma <- function(n, p=n, shape=1, scale=1){
	.CallOctave('randg', shape, n, p, scale)
}

#' Drawing from R Poisson Distribution in Octave
#' 
#' This function wraps a call to the standard Octave function \code{randp}, 
#' which is redefined by \code{RcppOctave} to call the R base function 
#' \code{\link[stats]{rpois}}.
#' This enables to exactly reproduce stochastic computations in R and Octave, 
#' without changing the original Octave/Matlab code.
#' See \code{\link{o_runif}} for more details.
#' 
#' This function calls the Octave function \code{randp} as provided by 
#' \code{RcppOctave}, which returns the draws as double.
#' The result is converted into integer to match the behaviour of
#' \code{\link{rpois}}
#' 
#' @inheritParams o_runif 
#' @param lambda Mean of the Poisson distribution
#' 
#' @templateVar name randp
#' @template OctaveDoc
#' 
#' @export
#' @family orandom
#' @seealso rpois 
#' @examples
#' 
#' # Draw random poisson values (in vector form)
#' set.seed(123)
#' o_rpois(1, lambda = 4)
#' o_rpois(1, 10, 4)
#' 
#' # Draw random poisson values (in matrix form)
#' set.seed(123)
#' o_rpois(2, lambda = 4)
#' o_rpois(2, 5, lambda = 4)
#' 
o_rpois <- function(n, p=n, lambda){
	res <- .CallOctave('randp', lambda, n, p)
    # need to convert to integer
    storage.mode(res) <- 'integer'
    res
}

##' Seed R RNG in Octave
##' 
##' This function seeds the Octave port of R random number generators.
##' 
##' The RNGs are provided by the Octave module \code{Rrng.oct}.
##' This module is built against the \code{libRrng} library, that provides R RNG 
##' original implementation, extracted from R-2.14.0. 
##' 
##' @param seed a single integer seed.
##' @return None
##' 
##' @export
##' @seealso set.seed
##' @family orandom
##' @examples
##' 
##' # Perform stochastic computations is R
##' set.seed(12345)
##' runif(5)
##' runif(5)
##' \dontshow{ rr <- runif(10) }
##' 
##' # Reproduce in Octave
##' o_set.seed(12345)
##' o_runif(1, 5)
##' o_runif(1, 5)
##' \dontshow{ stopifnot( identical(rr, o_runif(1, 10)) ) }
##'  
##' o_set.seed(12345)
##' o_runif(5)
##' 
#o_set.seed <- function(seed){
#	.CallOctave('setseed', seed)
#	invisible()
#}

##' Get RNG in Octave
##' 
##' @return None
##' 
##' @export
##' @rdname o_set.seed
##' @examples
##' 
##' # Set RNG seed in R
##' set.seed(12345)
##' seedR <- .Random.seed
##' 
##' # Set RNG seed in Octave
##' o_set.seed(12345)
##' seedO <- o_Random.seed()
##' 
##' identical(seedR, seedO)
##' \dontshow{ stopifnot( identical(seedR, seedO) ) }
##' 
##' identical( runif(10), o_runif(1, 10) )
##' \dontshow{ stopifnot( identical( runif(10), o_runif(1, 10) ) ) }
##'  
#o_Random.seed <- function(){
#	.CallOctave('getseed')	
#}


