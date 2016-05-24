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



# ++++++++++++++++++++++++++++++++++++++++++
# Functions shared by all the R2Cuba functions
# ++++++++++++++++++++++++++++++++++++++++++
verif <- function(ndim, ncomp, lower, upper, rel.tol, abs.tol,
            flags, min.eval,  max.eval) {
  # Verification of the input arguments
  # Issue a warning for each error and return T if none
  bon <- TRUE

  if ( (ndim >40) || (ndim<=0)) {
    bon <- FALSE
      warning("ndim should be positive and less or equal to 40")
  }
  if ( (ncomp >10) || (ncomp<=0)) {
    bon <- FALSE
      warning("ncomp should be positive and less or equal to 10")
  }

  if ((length(lower) != ndim) ||
      (length(upper) != ndim)) {
       bon <- FALSE
      warning("lower and upper should be vectors of length ndim")
  }

    if (any(lower >= upper)) {
      bon <- FALSE
      warning("Lower bounds should be less than upper bounds")
    }
  if ( !is.null(flags$verbose) &&
      ((flags$verbose <0) || (flags$verbose >3))) {
     bon <- FALSE
      warning("flags$verbose should be in [0,3]")
   }

  if ( !is.null(flags$final) &&
       ((flags$final <0) || (flags$final >1))) {
     bon <- FALSE
      warning("flags$final should be in [0,1]")
   }
  
      if ( !is.null(flags$pseudo.random) &&
       ((flags$pseudo.random <0) || (flags$pseudo.random >1))) {
     bon <- FALSE
      warning("flags$pseudo.random should be in [0,1]")
   }
      if ( !is.null(flags$smooth) &&
       ((flags$smooth <0) || (flags$smooth >1))) {
     bon <- FALSE
      warning("flags$smooth should be in [0,1]")
   }

  if ( (min.eval<0) || (min.eval>max.eval)) {
    bon <- FALSE
      warning("Error in min.eval or max.eval")
   }
  
     return(bon)
  } # End verif
# -------------------------------------------
decodflags <- function(flags) {
# Decode the flags
    if (!is.null(flags$verbose))
      lesflags <- flags$verbose
    else
      lesflags <- 1 # default value
    if (!is.null(flags$final))
      lesflags <- lesflags + (flags$final*4)
    if (!is.null(flags$pseudo.random))
      lesflags <- lesflags + (flags$pseudo.random*8)
    if (!is.null(flags$smooth))
      lesflags <- lesflags + (flags$smooth*16)
return(lesflags)
  } # End decodflags

# -------------------------------------------------------

crff <- function(lecall, integrand, nomf, libargs, ...) {
  # Determine how to call the user function according to
  # the list of its arguments and the current list of arguments

 #nfarg: number of the R user function formals arguments
 nfarg <- length(formals(integrand))
 
 if (nfarg <1)
   stop("Function integrand should have one argument at least")
 
 # nargsup:  number of additional arguments in the current call
 zl <- names(lecall)[-1]
 a <- zl %in% libargs # arguments only for the integration algorithm
 nargsup <- length(a[a==FALSE])
 
 if (nargsup >0) {
   # call with additional arguments
   if (nfarg == (nargsup+1))
     ffintegrand <- function(x, phw=0) integrand(x, ...)
   else
     if (nfarg >= (nargsup+2)) {
       ffintegrand <- function(x, phw=0) integrand(x,phw, ...)
     }
   else {
     stop(paste("Additional argument", zl[a==FALSE], "not expected in the integrand function\n"))
   }

 } # End (nargsup >0)
 else {
  if (nfarg == 1)
     ffintegrand <- function(x, phw=0) integrand(x)
   else {
      ffintegrand <- function(x, phw=0) integrand(x,phw)
   }
 }
 return(ffintegrand)
} # End crff

# ++++++++++++++++++++++++++++++++++++++++++
# The methods of the class "cuba"
# ++++++++++++++++++++++++++++++++++++++++++

print.cuba <- 
          function(x,  ...)
  {
                for (i in 1:length(x$value)) {
       cat("integral: ", format(x$value[i], ...),
           " (+-",
            format(x$abs.error[i], digits = 2), ")\n", sep = "")
            if (!is.null(x$nregions))
              cat("nregions: ", x$nregions, "; ",sep="")
     }
            cat("number of evaluations: ", x$neval)
             cat("; probability: ", format(x$prob, ...), "\n")

              if (x$message !="OK") 
         cat("failed with message ", sQuote(x$message), "\n")
                if ((x$ifail >1) && (x$method=="divonne")) {
                cat(" Divonne has estimated the number of points by which max.eval needs to be increased to reach the desired accuracy. It is", x$ifail,  "\n")
              }
                  
    invisible(x)

          }

