#
# file: estimate1d.R
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
#' @include rdist90ci_exact.R 
#' @include rdistq_fit.R
#' @include rtnorm90ci.R 
NULL
##############################################################################################
# estimate1d(distribution, lower, upper,...)
##############################################################################################
#' Create a 1-dimensional estimate object.
#' 
#' \code{estimate1d} creates an object of class \code{estimate1d}. The estimate of a one dimensional
#' variable is at minimum defined by the type of a univariate parametric distribution, the 5\% - and
#' 95\% quantiles. Optionally, the median can be supplied.
#' @param distribution \code{character}: A character string that defines the type of the univariate
#'   parametric distribution. 
#' @param lower \code{numeric}: lower bound of the 90\% confidence interval, i.e the 5\%-quantile 
#'   of this estimate.
#' @param upper \code{numeric}: upper bound of the 90\% confidence interval, i.e the 95\%-quantile 
#'   of this estimate.
#' @param ... arguments that can be coerced to a list comprising further elements of the 1-d 
#'   estimate (for details cf. below). Each element must be atomic and of length 1 (1-d property).
#' @details 
#'   It must hold that \code{lower <= upper}.
#'   \subsection{The structure of the input arguments}{
#'     \subsection{Mandatory input elements}{
#'     \tabular{lll}{
#'       \bold{Argument}    \tab  \bold{R-type}     \tab \bold{Explanation}\cr
#'       \code{distribution} \tab  \code{character} \tab  Distribution type of the estimate \cr
#'       \code{lower}        \tab  \code{numeric}   \tab   5\%-quantile of the estimate\cr
#'       \code{upper}        \tab  \code{numeric}   \tab  95\%-quantile of the estimate
#'     }
#'     }
#'     \subsection{Optional input elements}{
#'     The optional parameters in \code{...} provide additional characteristics of the 1-d estimate. 
#'     Frequent optional elements are:
#'     \tabular{lll}{
#'       \bold{Argument}     \tab  \bold{R-type}                 \tab \bold{Explanation}\cr
#'       \code{variable}     \tab  \code{character}              \tab  Variable name\cr
#'       \code{median}       \tab  cf. below                     \tab  50\%-quantile of the estimate\cr
#'       \code{method}       \tab  \code{character}              \tab  Method for calculation of distribution parameters
#'    }
#'    \subsection{The \code{median}}{
#'      If supplied as input, \code{median} can be either \code{NULL},  \code{numeric} or the
#'      character string \code{"mean"}. If it is \code{NA} it is set to \code{NULL}; if it equals
#'      \code{"mean"} it is set to \code{mean(c(lower, upper))}; if it is \code{numeric} it must 
#'      hold that \code{lower <= median <= upper}. 
#'      In case that no element \code{median} is provided, the default is \code{median=NULL}. 
#'    }    
#'    }
#'  }
#' @return An object of class \code{estimate1d} and \code{list} with at least (!) the elements:
#'    \tabular{lll}{
#'      \bold{Element}      \tab  \bold{R-type}                 \tab \bold{Explanation}\cr
#'      \code{distribution} \tab  \code{character}              \tab  Distribution type of the estimate \cr
#'      \code{lower}        \tab  \code{numeric}                \tab   5\%-quantile of the estimate\cr
#'      \code{median}       \tab  \code{numeric} or \code{NULL} \tab  50\%-quantile of the estimate\cr
#'      \code{upper}        \tab  \code{numeric}                \tab  95\%-quantile of the estimate
#'    }
#'    Note that the \emph{\code{median}} is a mandatory element of an \code{estimate1d}, although it
#'     is not necessary as input. If \code{median} is numeric it holds that:
#'    \code{lower <= median <= upper}. In any case an \code{estimate1d} object has the property 
#'    \code{lower <= upper}.
#' @seealso  \code{\link{random.estimate1d}}
#' @export
estimate1d<-function(distribution, lower, upper, ...){
  # Check preconditions:
  ## Mandatory arguments:
  ### Check argument types:
  if ( is.null(distribution) || !is.character(distribution))
    stop("\"distribution\" must be supplied as character string.")
  if ( is.null(lower) || is.na(lower<-as.numeric(lower)) )
    stop("\"lower\" must be supplied as numeric.")
  if ( is.null(upper) || is.na(upper<-as.numeric(upper)) )
    stop("\"upper\" must be supplied as numeric.")
  #### Check 1-d property:
  if ( length(distribution) > 1 || length(lower) > 1 || length(upper) > 1)
    stop("Input dimension is larger than 1. All input elements must be 1 dimensional.")
  ### Check input semantics:
  if ( lower > upper )
    stop("\"lower > upper\"")
  ## Optional arguments:
  optionalArguments<-if( missing(...) || length(unlist(list(...)))==0 ) NULL else as.list(unlist(list(...))) 
  #optionalArguments<-list(...)
  median<-NULL
  if( !is.null(optionalArguments) )
    for ( i in names(optionalArguments) ){
      ### Check argument types:
      if (i == "variable" && !is.character(optionalArguments[[i]]))
        stop("Optional argument \"", i, "\" is not character.")
      if (i == "method" && !is.character(optionalArguments[[i]]))
        stop("Optional argument \"", i, "\" is not character.")
      #Process median:
      if (i == "median"){
        median<-optionalArguments[[i]]
        optionalArguments<-optionalArguments[!names(optionalArguments) %in% "median"]
        if ( !is.null(median) ){
          if( is.na(median) )
            median<-NULL
          else if( is.character(median) && median=="mean")
            median<-as.numeric(mean(c(lower, upper)))
          #else if ( !is.numeric(median) ) stop("\"median=", median, "\" is not allowed!")
          else if (  is.na(median<-as.numeric(median)) ) 
            stop("\"median=", median, "\" is not allowed!")
          else if(  (lower > median || median > upper) )
            stop("It must hold: \"lower <= median <= upper\"")
        }
      }
      ### Check 1-d property:
      if ( !is.atomic(optionalArguments[[i]]) )
        stop("Optional argument \"", i, "\" is not atomic.")
      if ( length(optionalArguments[[i]]) > 1 )
        stop("Optional argument \"", i, "\" is not 1 dimensional.")
    }
  # Create the estimate1d:
  if( as.logical(length(optionalArguments)) )
    estimate1dObject<-as.list( c(distribution=distribution, lower=lower, median=median, upper=upper, optionalArguments) )
  else
    estimate1dObject<-list(distribution=distribution, lower=lower, median=median, upper=upper )    
  # Return object:
  class(estimate1dObject)<-c("estimate1d", class(estimate1dObject))
  estimate1dObject
}
##############################################################################################
# as.estimate1d(x, ...)
##############################################################################################
#' Transform to an 1-dimensional estimate object
#'
#' \code{as.estimate1d} tries to transform an object to class \code{estimate1d}.
#' @param x an object to be transformed to class \code{estimate1d}. 
#' @rdname estimate1d
#' @export
as.estimate1d<-function(x, ...){
  x_vec<-unlist(x)
  if ( !"distribution" %in% names(x_vec) )
    stop( "no distribution element!")
  if ( !"lower" %in% names(x_vec) )
    stop( "no lower element!")
  if ( !"upper" %in% names(x_vec) )
    stop( "no upper element!")
  estimate1dObject<-estimate1d(distribution=x_vec[["distribution"]], 
                               lower=x_vec[["lower"]],
                               upper=x_vec[["upper"]],
                               as.list(x_vec[!names(x_vec) %in% c("distribution", "lower", "upper")]))
  # Return object:
  estimate1dObject
}
##############################################################################################
# random.estimate1d(rho,n,method, ...)
##############################################################################################
#' Generate univariate random numbers defined by a 1-d estimate.
#' 
#' This function generates random numbers for univariate parametric distributions, whose  
#' parameters are determined by a one dimensional estimate (\code{\link{estimate1d}}).
#' @param rho \code{estimate1d}: Univariate distribution to be randomly sampled. 
#' @param n \code{integer}: Number of observations to be generated
#' @param method \code{character}: Particular method to be used for random number generation. It 
#'    can be either \code{"calculate"} (the default) or \code{"fit"}. Details below.
#' @param relativeTolerance \code{numeric}: the relative tolerance level of deviation of the
#'   generated confidence interval from the specified interval. If this deviation is greater than
#'   \code{relativeTolerance} a warning is given.
#' @param ... Optional arguments to be passed to the particular random number
#'  generating function (cf. below).
#'  @details
#'  \describe{
#'    \item{\code{rho[["distribution"]]}:}{
#'    The follwing table shows the available distributions and the implemented generation method:
#'    \tabular{lll}{
#'    \bold{\code{rho[["distribution"]]}}  \tab\bold{Distribution Name}                         \tab \bold{\code{method}} \cr
#'    \code{"const"}                \tab Deterministic case                            \tab not applicable\cr
#'    \code{"norm"}                 \tab \link{Normal}                                 \tab \code{\link[=rdist90ci_exact]{calculate}}, \code{\link[=rdistq_fit]{fit}}  \cr
#'    \code{"posnorm"}              \tab \link[=rposnorm90ci]{Positive normal}         \tab \code{\link[=paramtnormci_numeric]{calculate}}, \code{\link[=paramtnormci_fit]{fit}} \cr
#'    \code{"tnorm_0_1"}            \tab \link[=rtnorm_0_1_90ci]{0-1-truncated normal} \tab \code{\link[=paramtnormci_numeric]{calculate}}, \code{\link[=paramtnormci_fit]{fit}} \cr
#'    \code{"beta"}                 \tab \link{Beta}                                   \tab \code{\link[=rdistq_fit]{fit}}  \cr
#'    \code{"cauchy"}               \tab \link{Cauchy}                                 \tab \code{\link[=rdistq_fit]{fit}}  \cr
#'    \code{"logis"}                \tab \link{Logistic}                               \tab \code{\link[=rdistq_fit]{fit}}  \cr
#'    \code{"t"}                    \tab \link[=TDist]{Student t}                      \tab \code{\link[=rdistq_fit]{fit}}  \cr
#'    \code{"chisq"}                \tab \link[=Chisquare]{Central Chi-Squared}        \tab \code{\link[=rdistq_fit]{fit}}  \cr
#'    \code{"chisqnc"}              \tab \link[=Chisquare]{Non-central Chi-Squared}    \tab \code{\link[=rdistq_fit]{fit}}  \cr
#'    \code{"exp"}                  \tab \link{Exponential}                            \tab \code{\link[=rdistq_fit]{fit}}  \cr  
#'    \code{"f"}                    \tab \link[=FDist]{Central F}                      \tab \code{\link[=rdistq_fit]{fit}}  \cr
#'    \code{"gamma"}                \tab \link[=GammaDist]{Gamma} with \code{scale=1/rate} \tab \code{\link[=rdistq_fit]{fit}}  \cr
#'    \code{"lnorm"}                \tab \link[=Lognormal]{Log Normal}                 \tab \code{\link[=rdist90ci_exact]{calculate}}, \code{\link[=rdistq_fit]{fit}}  \cr
#'    \code{"unif"}                 \tab \link{Uniform}                                \tab \code{\link[=rdist90ci_exact]{calculate}}, \code{\link[=rdistq_fit]{fit}}  \cr
#'    \code{"weibull"}              \tab \link{Weibull}                                \tab \code{\link[=rdistq_fit]{fit}}  \cr
#'    \code{"triang"}               \tab \link[mc2d:triangular]{Triangular}            \tab \code{\link[=rdistq_fit]{fit}}  \cr
#'    \code{"gompertz"}             \tab \link[eha:Gompertz]{Gompertz}                 \tab \code{\link[=rdistq_fit]{fit}}  \cr
#'    \code{"pert"}                 \tab  \link[mc2d:pert]{(Modified) PERT}            \tab \code{\link[=rdistq_fit]{fit}}  \cr
#    \code{"tnorm"}                \tab  \link[msm:tnorm]{Truncated Normal}           \tab \code{\link[=rdistq_fit]{fit}} 
#'    }
#'    For \code{distribution="const"} the argument \code{method} is obsolete, as a constant is neither
#'    fitted nor calculated.  
#'    }
#'    \item{\code{rho[["method"]]}}{
#'    If supplied, i.e. \code{!is.null(rho[["method"]])}, this value overwrites the function 
#'    argument \code{method}. 
#'    }
#'    \item{\code{method}}{
#'    This parameter defines, how the parameters of the distribution to be sampled are derived from 
#'    \code{rho[["lower"]]}, \code{rho[["upper"]]} and possibly \code{rho[["median"]]}.
#'    Possibilities are \code{"calculate"} (the default) or \code{"fit"}:
#'    \describe{
#'      \item{\code{method="calculate"}}{
#'      The parameters are calculated if possible using the exact (analytical) formula or, otherwise,
#'      numerically. This calculation of the distribution parameters is independent of 
#'      \code{rho[["median"]]} being  supplied or not. For the implemented distributions, it only 
#'      depends on \code{rho[["lower"]]} and \code{rho[["upper"]]}. However, if it is supplied, i.e. 
#'      \code{is.numeric(rho[["median"]])}, a check is performed, if the relative deviation of the 
#'      generated median from \code{rho[["median"]]} is greater than \code{relativeTolerance}. In 
#'      this case a warning is given.
#'      }
#'      \item{\code{method="fit"}}{
#'      The parameters are obtained by fitting the distribution on the supplied quantiles.
#'      Given that \code{rho[["median"]]==NULL} the distribution is fitted only to \code{lower} and 
#'      \code{upper} and a warning is given; due to the used numerical procedure, the calculated 
#'      parameters might define a distribution which strongly deviates from the intended one. There is 
#'      larger control on the shape of the distribution to be generated by supplying the estimate of the 
#'      median. If \code{is.numeric(rho[["median"]])} the distribution is fitted to \code{lower}, 
#'      \code{upper} and \code{median}. 
#'      }   
#'    }
#'    }
#'    \item{\code{...}}{
#'      For passing further parameters to the function which generates the random numbers, cf. 
#'      the above table and follow the link in the column \code{method}.
#'    }
#'  } 
#'  
#' @seealso \code{\link{estimate1d}}; For \code{method="calculate"}: \code{\link{rdist90ci_exact}}; for \code{method="fit"}: \code{\link{rdistq_fit}}; for both
#'   methods: \code{\link{rposnorm90ci}} and \code{\link{rtnorm_0_1_90ci}}. For the default method: \code{\link{random}}.

#' @export
random.estimate1d<-function(rho ,n , method="calculate", relativeTolerance=0.05, ...){
  # Overwrite argument "method" if rho[["method"]] is supplied:
  if ( !is.null(rho[["method"]]) && rho[["method"]] !="")
    method<-rho[["method"]]
  # Create output vector for the random numbers to be generated:
  x<-vector(length=n)
  # Generate the random numbers according to the distribution type:
  ## Constants are neither calculated nor fitted, i.e. the procedure is the same for all methods as they are constant:
  if(match(rho[["distribution"]], "const", nomatch = 0)){
    x <-  rdist90ci_exact(distribution="const",
                          n=n,
                          lower=rho[["lower"]],
                          upper=rho[["upper"]])
  } 
  ## Generate the random numbers by calculating the distribution parameters from lower and upper:
  else if(method=="calculate"){
    # ToDo: extract this block as function rdist90ci_calculate()?
    if(match(rho[["distribution"]], c("norm", 
                                      "lnorm",
                                      "unif"), nomatch = 0)){
      x <-  rdist90ci_exact(distribution=rho[["distribution"]],
                            n=n,
                            lower=rho[["lower"]],
                            upper=rho[["upper"]])
    }
    else if(match(rho[["distribution"]], c("posnorm"), nomatch = 0)){
      x <-  rposnorm90ci(n=n,
                         lower=rho[["lower"]],
                         upper=rho[["upper"]],
                         method="numeric",
                         relativeTolerance = relativeTolerance)
    } 
    else if(match(rho[["distribution"]], c("tnorm_0_1"), nomatch = 0)){
      x <-  rtnorm_0_1_90ci(n=n,
                            lower=rho[["lower"]],
                            upper=rho[["upper"]],
                            method="numeric",
                            relativeTolerance = relativeTolerance)
    }
    else if(match(rho[["distribution"]], c("triang"), nomatch = 0)){
      if(is.na(rho[["median"]]))
        stop ("Median missing for calculation of", rho[["distribution"]], "distribution.")
      x <-  mc2d::rtriang(n=n,
                    min=rho[["lower"]],
                    mode=rho[["median"]],
                    max=rho[["upper"]])
      warning("For the triangular distribution, the lower and upper bounds are interpreted as the extreme values (not the 90% confidence interval) and the 'median' as the mode (the most likely value).")
      
    }
    
    else
      stop("\"", rho[["distribution"]], "\" is not a valid distribution type for method=\"", method, "\".")
    ### Generate warning if relative deviation from median (if provided) is greater than 
    ### relativeTolerance:
    if(!rho[["distribution"]]=="triang")
      if ( !is.null(median<-rho[["median"]]) ){
        median_calc<- quantile(x=x,probs=0.50)[["50%"]]
        scale <- if( median > 0 ) median else NULL
        if( !isTRUE( msg<-all.equal(median, median_calc,  scale=scale, tolerance=relativeTolerance) ) ){
          warning("For method=\"calculate\": deviation of calculated \"median\" from supplied target value:\n",
                  "  Calculated value: ", median_calc, "\n  ",
                  "Target value:     ", median,   "\n  ",
                  msg)
        }
      }
  }
  ## Generate the random numbers by fitting the distribution parameters on lower and upper, and 
  ##  (if provided) on the median:
  else if (method=="fit"){
    if ( is.null(rho[["median"]]) || rho[["distribution"]]=="unif"){
      percentiles<-c(0.05,0.95)
      quantiles<-c(rho[["lower"]], rho[["upper"]])
    } else {
      percentiles<-c(0.05,0.5,0.95)
      quantiles<-c(rho[["lower"]], rho[["median"]], rho[["upper"]])
    }
    if(match(rho[["distribution"]], c("norm", 
                                      "beta",
                                      "cauchy",
                                      "logis",
                                      "t",
                                      "chisq",
                                      "chisqnc",
                                      "exp",
                                      "f",
                                      "gamma",
                                      "lnorm",
                                      "unif",    
                                      "weibull",
                                      "triang",
                                      "gompertz"), nomatch = 0)){
      x<-rdistq_fit(distribution=rho[["distribution"]], 
                    n=n, 
                    percentiles=percentiles, 
                    quantiles=as.numeric(quantiles), 
                    relativeTolerance=relativeTolerance,
                    ...) 
    }  
    else if(match(rho[["distribution"]], c("posnorm"), nomatch = 0)){
      x <-  rposnorm90ci(n=n,
                         lower=rho[["lower"]],
                         median=rho[["median"]],
                         upper=rho[["upper"]],
                         method="fit",
                         relativeTolerance=relativeTolerance)
    } 
    else if(match(rho[["distribution"]], c("tnorm_0_1"), nomatch = 0)){
      x <-  rtnorm_0_1_90ci(n=n,
                            lower=rho[["lower"]],
                            median=rho[["median"]],
                            upper=rho[["upper"]],
                            method="fit",
                            relativeTolerance=relativeTolerance)
    } 
    else
      stop("\"", rho[["distribution"]], "\" is not a valid distribution type for method=\"", method, "\".")
    ### Generate warning if median is not provided:
    if ( is.null(rho[["median"]]) )
      warning("For method=\"fit\": no \"median\" supplied. Sometimes, this might not lead to the the desired
              distributions shape. Advise: check it.")
  }
  else
    stop ("method must be either \"calculate\" or \"fit\".")
  # Return generated random numbers:
  x
}

#' @examples
#' # Generate log normal distributed random numbers:
#' x<-random(estimate1d("lnorm",50,100), n=100000)
#' quantile(x, probs=c(0.05, 0.95))
#' hist(x, breaks=100)
