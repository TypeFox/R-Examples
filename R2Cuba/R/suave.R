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

suave <- function(ndim, ncomp,
                  integrand, ...,
                    lower=rep(0,ndim), upper=rep(1,ndim),
                   rel.tol= 1e-3,
                  abs.tol = 0,
                  flags=list(verbose=1,
                  final=1, pseudo.random=0, smooth=0,
                    mersenne.seed=NULL),
                  min.eval=0,  max.eval=50000,
                  nnew=1000, flatness= 50)
  {
     # Verification
 if (!verif(ndim, ncomp, lower, upper, rel.tol, abs.tol,
            flags, min.eval,  max.eval))
      stop("Error in input: see the warnings")

 if (nnew <= 0)
   stop("nnew should be positive")
 
  # Decode the flags
    lesflags <- decodflags(flags)
 if (is.null(flags$mersenne.seed))
  flags$mersenne.seed <- NA
                                           
    # Allocate output structures
    nregions <- 0
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
               "nnew", "flatness")
 ffintegrand <- crff(match.call(), integrand, "suave", libargs, ...)

 prdbounds <- prod(upper-lower)
    
  ret <-  .C("Rsuave", as.integer(ndim),
             as.integer(ncomp),
             ffintegrand, new.env(),
             as.double(lower), as.double(upper),
             as.double(prdbounds),
             as.double(rel.tol), as.double(abs.tol),
              as.integer(flags$mersenne.seed),
             as.integer(lesflags),
             as.integer(min.eval),  as.integer(max.eval),
             as.integer(nnew),  as.double(flatness),
             nregions=as.integer(nregions),
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
      mess ="OK" # OK required to be printed by print.cuba
    else
      mess="Dimension out of range"

    retour <- list(method="suave",
                   neval=ret$neval, nregions=ret$nregions,
                   ifail=ret$fail, value=ret$integral,
                abs.error=ret$error,
                   neval=ret$neval,
                   prob=ret$prob, message=mess,
                call=match.call())
    attr(retour, "class") = c("cuba")
    return(retour)
  } # End suave

                
