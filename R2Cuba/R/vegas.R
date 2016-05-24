###################################################################
# R2Cuba R package
# Copyright INRA 2015
# INRA, UR1404, Research Unit MaIAGE
# F78352 Jouy-en-Josas, France.
#
# This file is part of R2Cuba R package, interface between R and
# the Cuba library.
# R2Cuba URL: http://cran.r-project.org/web/packages/R2Cuba
# Cuba library URL: http://www.feynarts.de/cuba

# R2Cuba is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See the GNU General Public License at:
# http://www.gnu.org/licenses/
#
###################################################################

# -------------------------------------------  
verifvegas <- function( nstart, nincrease, nbatch,gridno, state.file ) {
  # Verification
  bon <- TRUE
  if ( !is.null(nbatch) &&
      (nbatch<=0)) {
     bon <- FALSE
      warning("nbatch should be positive")
   }
  if ( !is.null(gridno) &&
      (gridno > 10)) {
     bon <- FALSE
      warning("gridno less or equal to 10")
   }
  
if (nstart <=0) {
  bon <- FALSE
      warning("nstart should be positive")
}

if (nincrease <=0) {
  bon <- FALSE
      warning("nincrease should be positive")
}

  if (!is.null(state.file)) {
    if ( !is.character(state.file) ||  (nchar(state.file)>= 128)) {
    bon <- FALSE
      warning("state.file should be a character string of length less than 128")
  }
  }
  
    
  return(bon)
  } # End verifvegas



# ____________________________________
vegas <- function(ndim, ncomp,
                integrand, ...,
                 lower=rep(0,ndim), upper=rep(1,ndim),
                  rel.tol= 1e-3,
                  abs.tol = 0,
                   flags=list(verbose=1,
                  final=1, pseudo.random=0, smooth=0,
                     mersenne.seed=NULL),
                  min.eval=0,  max.eval=50000,
                  nstart=1000,  nincrease=500,
                  nbatch=1000, gridno=0, state.file=NULL)
  {
     # Verification
if (!verif(ndim, ncomp, lower, upper, rel.tol, abs.tol,
            flags, min.eval,  max.eval) ||
    !verifvegas( nstart, nincrease, nbatch, gridno, state.file))
      stop("Error in input: see the warnings")


 # Decode the flags
    lesflags <- decodflags(flags)
if (is.null(flags$mersenne.seed))
  flags$mersenne.seed <- NA
if (is.null(nbatch))
  nbatch <- 1000
if (is.null(gridno))
  gridno <- 0
if (is.null(state.file))
  state.file <- ""

    # Allocate output structures
    neval <- 0
    fail <- 0
    integral <- rep(0, ncomp)
    error <- rep(0, ncomp)
    prob <- rep(0, ncomp)
  #  ffintegrand <-
#	if(length(list(...)) && length(formals(integrand)) > 2)
#	    function(x, weight) integrand(x, weight, ...)
#	else integrand
  libargs <- c("ndim", "ncomp", 
                 "integrand","lower", "upper",
                  "rel.tol", "abs.tol", "flags",
                    "min.eval","max.eval", 
                  "nstart",  "nincrease",
                  "nbatch", "gridno", "state.file")
 ffintegrand <- crff(match.call(), integrand, "vegas", libargs, ...)

    
  prdbounds <- prod(upper-lower)
  ret <-  .C("Rvegas", as.integer(ndim),
             as.integer(ncomp),
             ffintegrand, new.env(),
             as.double(lower), as.double(upper),
             as.double(prdbounds),
             as.double(rel.tol), as.double(abs.tol),
             as.integer(flags$mersenne.seed),
             as.integer(nbatch),
             as.integer(gridno), 
             as.integer(lesflags),
             as.integer(min.eval),  as.integer(max.eval),
             as.integer(nstart),  as.integer(nincrease),
             as.character(state.file),
             neval=as.integer(neval), fail=as.integer(fail),
             integral=as.double(integral),
             error=as.double(error), prob=as.double(prob),
             NAOK=TRUE)
#Add to finish the last print:cat("\n")

# To homogeneize with the R function "integrate", add
    # message and call into the output,
    # ifail rather than fail , abs.error rather than  error,
    # value rather than integral
    
    
    if (ret$fail ==0)
      mess ="OK" # compulsory for the printing by print.cuba
    else
      mess="Dimension out of range"

    retour <- list(method="vegas",
                   neval=ret$neval,
                   ifail=ret$fail, value=ret$integral,
                abs.error=ret$error,
                   neval=ret$neval,
                   prob=ret$prob, message=mess,
                call=match.call())
    attr(retour, "class") = c("cuba")
    
 
    return(retour)
  } # End vegas

                
