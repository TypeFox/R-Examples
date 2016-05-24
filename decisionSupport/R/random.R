# 
# file: random.R
#
# This file is part of the R-package decisionSupport
# 
# Authors: 
#   Lutz Göhring <lutz.goehring@gmx.de>
#   Eike Luedeling (ICRAF) <eike@eikeluedeling.com>
#
# Copyright (C) 2015 World Agroforestry Centre (ICRAF) 
#	http://www.worldagroforestry.org
# 
# The R-package decisionSupport is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# The R-package decisionSupport is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with the R-package decisionSupport.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################################
#' @include rdistq_fit.R
NULL
##############################################################################################
# generic: random(rho,n,method,...)
# ToD: rename to rdistq or rq?
##############################################################################################
#' Quantiles or empirically based generic random number generation.
#' 
#' These functions generate random numbers for parametric distributions, parameters of which are 
#' determined by given quantiles or for distributions purely defined empirically. 
#' @param rho Distribution to be randomly sampled.
#' @param n \code{integer}: Number of observations to be generated
#' @param method \code{character}: Particular method to be used for random number generation. 
#' @param relativeTolerance \code{numeric}: the relative tolerance level of deviation of the
#'   generated confidence interval from the specified interval. If this deviation is greater than
#'   \code{relativeTolerance} a warning is given.
#' @param ... Optional arguments to be passed to the particular random number
#'  generating function.
#' @export
random <- function(rho,n,method, relativeTolerance, ...) UseMethod("random")
##############################################################################################
# random.default(rho,n,method,...)
##############################################################################################
#' The default method generates univariate random numbers specified by arbitrary quantiles.
#' @describeIn random Quantiles based univariate random number generation.
#'  \describe{
#'    \item{\bold{Arguments}}{
#'       \describe{  
#'       \item{\code{rho} }{
#'         rho \code{list}: Distribution to be randomly sampled. The list elements are 
#'          \code{$distribution}, \code{$probabilities} and \code{$quantiles}. For details cf. below.
#'       }
#'       \item{\code{method} }{
#'         \code{character}: Particular method to be used for random number 
#'         generation. Currently only method \code{\link{rdistq_fit}{fit}} is implemented which is the 
#'         default.
#'       }
#'       \item{\code{relativeTolerance}}{
#'         \code{numeric}: the relative tolerance level of deviation of the generated confidence 
#'         interval from the specified interval. If this deviation is greater than 
#'         \code{relativeTolerance} a warning is given.
#'       }
#'        \item{\code{...}}{
#'          Optional arguments to be passed to the particular random number
#'          generating function, i.e. \code{\link{rdistq_fit}}.
#'        }
#'      }
#'    }
#'    \item{\bold{Details}}{
#'      \describe{
#'        \item{ }{
#'          The distribution family is determined by \code{rho[["distribution"]]}. For the  
#'          possibilities cf. \code{\link{rdistq_fit}}.
#'        }
#'        \item{ }{
#'          \code{rho[["probabilities"]]} and \code{[[rho"quantiles"]]} are numeric vectors of the same 
#'          length. The first defines the probabilities of the quantiles, the second defines the quantiles 
#'          values which determine the parametric distribution.
#'        }
#'      }
#'    }
#'    \item{\bold{Value}}{
#'      \describe{
#'        \item{ }{
#'          A numeric vector of length \code{n} containing the generated random numbers.
#'        }
#'      }
#'    }
#'    \item{\bold{See Also}}{
#'      \describe{
#'        \item{ }{
#'          \code{\link{rdistq_fit}}
#'        }  
#'      } 
#'    } 
#'  }    
#' @examples
#'  x<-random(n=10000)
#'  hist(x,breaks=100)
#'  mean(x)
#'  sd(x)
#'   
#'  rho<-list(distribution="norm", 
#'            probabilities=c(0.05,0.4,0.8), 
#'            quantiles=c(-4, 20, 100))
#'  x<-random(rho=rho, n=10000, tolConv=0.01)
#'  hist(x,breaks=100)
#'  quantile(x,p=rho[["probabilities"]])
#' @export
random.default <- function(rho=list(distribution="norm", probabilities=c(0.05,0.95), quantiles=c( -qnorm(0.95), qnorm(0.95) )), 
                           n, method="fit", relativeTolerance=0.05, ...){
  # Check preconditions:
  if ( !is.list(rho) )
    stop("rho must be a list with elements \"distribution\", \"probabilities\" and \"quantiles\"")
  if( is.null(rho[["distribution"]]) )
    stop("rho[[\"distribution\"]] must be supplied.")
  if( is.null(rho[["probabilities"]]) || !all(!is.na(as.numeric(rho[["probabilities"]]))) )
    stop("rho[\"probabilities\"] must be supplied.")
  if( is.null(rho[["quantiles"]]) || !all(!is.na(as.numeric(rho[["quantiles"]]))) )
    stop("rho[\"quantiles\"] must be supplied.")
  if( length(rho[["probabilities"]])!=length(rho[["quantiles"]]) )
    stop( "length(rho[[\"probabilities\"]])!=length(rho[[\"quantiles\"]])" )
  # Constants are neither calculated nor fitted, i.e. the procedure is the same for all methods as they are constant:
  if(match(rho["distribution"], "const", nomatch = 0)){
    stop("const not implemented, yet")
  } 
  else if (method=="fit"){
    # The next few lines apply a curve fitting procedure based on given distributions and specified quantiles:
    if(match(rho["distribution"], c("norm", 
                                    "beta",
                                    "cauchy",
                                    "logis",
                                    "t",
                                    "chisq",
                                    "exp",
                                    "f",
                                    "gamma",
                                    "lnorm",
                                    "unif",    
                                    "weibull",
                                    "triang",
                                    "gompertz"), nomatch = 0)){
      x<-rdistq_fit(distribution=rho["distribution"], 
                    n=n, 
                    percentiles=as.numeric(rho[["probabilities"]]), 
                    quantiles=as.numeric(rho[["quantiles"]]), 
                    relativeTolerance=relativeTolerance,
                    ...) 
    }
    else
      stop("\"", rho[["distribution"]], "\" is not a valid distribution type for method=\"", method, "\".")
  }
  else
    stop ("method must be \"fit\".")
  # Return generated random numbers:
  x
}

##############################################################################################
# random.vector(rho,n,method,...)
##############################################################################################
#' \code{random.vector} generates univariate random numbers drawn from a distribution purely defined
#' empirically. 
#' @describeIn random Univariate random number generation by drawing from a given 
#'    empirical sample.
#'    \describe{
#'       \item{\bold{Arguments}}{
#'         \describe{
#'           \item{\code{rho} }{
#'             \code{vector}: Univariate empirical sample to be sampled from.
#'            }
#'            \item{\code{method} }{ 
#'               for this class no impact
#'            }
#'            \item{\code{relativeTolerance}}{
#'              for this class no impact
#'             }
#'            \item{\code{...}}{
#'              for this class no impact
#'            }
#'         }
#'      }
#'      \item{\bold{Value}}{
#'        \describe{
#'          \item{ }{
#'             A \code{numeric vector} of length \code{n} containing the generated random numbers.
#'          }
#'        }
#'      }
#'      \item{\bold{See Also}}{
#'        \describe{
#'           \item{ }{
#'             \code{\link{sample}}
#'           } 
#'        } 
#'      } 
#'    }  
#' @export
random.vector <- function(rho=runif(n=n), 
                          n, method=NULL, relativeTolerance=NULL, ...){
  # Return generated random numbers:
  sample(rho, size=n, replace=TRUE)
}
##############################################################################################
# random.data.frame(rho,n,method,...)
##############################################################################################
#' \code{random.data.frame} generates multivariate random numbers drawn from a distribution 
#' purely defined empirically. 
#' @describeIn random Multivariate random number generation by drawing from a given empirical sample.
#'   \describe{
#'     \item{\bold{Arguments}}{
#'       \describe{
#'         \item{\code{rho} }{
#'           \code{data.frame}: Multivariate empirical sample to be sampled from.
#'         }
#'         \item{\code{method} }{ 
#'           for this class no impact
#'         }
#'         \item{\code{relativeTolerance}}{
#'           for this class no impact
#'         }
#'          \item{\code{...}}{
#'            for this class no impact
#'         }
#'       }
#'     }
#'     \item{\bold{Value}}{
#'       \describe{
#'         \item{ }{
#'            A \code{data.frame} with \code{n} rows containing the generated random numbers.  
#'          }
#'        }
#'      }
#'      \item{\bold{See Also}}{
#'        \describe{
#'          \item{ }{
#'            \code{\link{sample}}
#'          }
#'        }  
#'      }
#'    }
#' @export
random.data.frame <- function(rho=data.frame(uniform=runif(n=n)), 
                              n, method=NULL, relativeTolerance=NULL, ...){
  # Return generated random numbers:
  data.frame(rho[sample.int(n=nrow(rho), size=n, replace=TRUE),], row.names=NULL)
}
