#
# file: rdistq_fit.R
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
##############################################################################################
# rdistq_fit(distribution, n, percentiles, quantiles)
##############################################################################################
#' Quantiles based univariate random number generation (by parameter fitting).
#' 
#' This function generates random numbers for a set of univariate parametric distributions from  
#' given quantiles. Internally, this is achieved by fitting the distribution function
#' to the given quantiles. 
#' @param distribution A character string that defines the univariate distribution
#'  to be randomly sampled. 
#' @param n Number of generated observations.
#' @param percentiles Numeric vector giving the percentiles. 
#' @param quantiles Numeric vector giving the quantiles. 
#' @param relativeTolerance \code{numeric}; the relative tolerance level of deviation of the
#'   generated individual percentiles from the specified percentiles. If any deviation is greater 
#'   than \code{relativeTolerance} a warning is given.
#' @inheritParams rriskDistributions::rriskFitdist.perc
#' @param verbosity \code{integer}; if \code{0} the function is silent; the larger the value the
#'   more verbose is the output information.
#' @details
#'   The following table shows the available distributions and their identification 
#'  (option: \code{distribution}) as a character string:
#'  \tabular{llll}{
#'  \bold{\code{distribution}}  \tab  \bold{Distribution Name}                 \tab \bold{\code{length(quantiles)}} \tab \bold{Necessary Package}\cr
#'  \code{"norm"}       \tab  \link{Normal}                                    \tab >=2    \tab \cr
#'  \code{"beta"}       \tab  \link{Beta}                                      \tab >=2    \tab \cr  
#'  \code{"cauchy"}     \tab  \link{Cauchy}                                    \tab >=2    \tab \cr
#'  \code{"logis"}      \tab  \link{Logistic}                                  \tab >=2    \tab \cr 
#'  \code{"t"}          \tab  \link[=TDist]{Student t}                         \tab >=1    \tab \cr 
#'  \code{"chisq"}      \tab  \link[=Chisquare]{Central Chi-Squared}           \tab >=1    \tab \cr 
#'  \code{"chisqnc"}    \tab  \link[=Chisquare]{Non-central Chi-Squared}       \tab >=2    \tab \cr 
#'  \code{"exp"}        \tab  \link{Exponential}                               \tab >=1    \tab \cr    
#'  \code{"f"}          \tab  \link[=FDist]{Central F}                         \tab >=2    \tab \cr 
#'  \code{"gamma"}      \tab  \link[=GammaDist]{Gamma} with \code{scale=1/rate}\tab >=2    \tab \cr
#'  \code{"lnorm"}      \tab  \link[=Lognormal]{Log Normal}                    \tab >=2    \tab \cr
#'  \code{"unif"}       \tab  \link{Uniform}                                   \tab ==2    \tab \cr
#'  \code{"weibull"}    \tab  \link{Weibull}                                   \tab >=2    \tab \cr
#'  \code{"triang"}     \tab  \link[mc2d:triangular]{Triangular}               \tab >=3    \tab \pkg{\code{mc2d}}\cr
#'  \code{"gompertz"}   \tab  \link[eha:Gompertz]{Gompertz}                    \tab >=2    \tab \pkg{\code{eha}} \cr
#'  \code{"pert"}       \tab  \link[mc2d:pert]{(Modified) PERT}                \tab >=4    \tab \pkg{\code{mc2d}}\cr
#'  \code{"tnorm"}      \tab  \link[msm:tnorm]{Truncated Normal}               \tab >=4    \tab \pkg{\code{msm}}
#'  }
#'  
#'  \code{percentiles} and \code{quantiles} must be of the same length. \code{percentiles} must be 
#'  \code{>=0} and \code{<=1}.
#'  
#'  
#' The default for \code{percentiles} is 0.05, 0.5 and 0.95, so for the default, 
#' the quantiles argument should be a vector with 3 elements. If this is to be longer,
#' the percentiles argument has to be adjusted to match the length of quantiles.
#' 
#'  The fitting of the distribution parameters is done using 
#'  \code{\link[rriskDistributions]{rriskFitdist.perc}}.
#' @return A numeric vector of length \code{n} with the sampled values according to the chosen 
#'   distribution.
#' @seealso \code{\link[rriskDistributions]{rriskFitdist.perc}}
#' @examples
#' # Fit a log normal distribution to 3 quantiles:
#' if ( requireNamespace("rriskDistributions", quietly = TRUE) ){
#'   percentiles<-c(0.05, 0.5, 0.95)
#'   quantiles=c(1,3,15)
#'   hist(r<-rdistq_fit(distribution="lnorm", n=10000, quantiles=quantiles),breaks=100)
#'   print(quantile(x=r, probs=percentiles))
#' }
#' @export
rdistq_fit <- function(distribution, n, percentiles=c(0.05,0.5,0.95), quantiles, 
                       relativeTolerance=0.05, tolConv=0.001, fit.weights=rep(1,length(percentiles)),
                       verbosity=1){
  # Check preconditions:
  ## Namespace requirements:
  if (!requireNamespace("rriskDistributions", quietly = TRUE)) 
    stop("Package \"rriskDistributions\" needed. Please install it, e.g. from 
         http://cran.r-project.org/src/contrib/Archive/rriskDistributions/",
         call. = FALSE)
  if(distribution=="triang"){
    requiredPackage<-"mc2d"
    if( !requireNamespace(requiredPackage, quietly = TRUE) ) 
      stop("Package \"",requiredPackage,"\" needed for option distribution=", distribution, ". Please install it.",
           call. = FALSE)
    else
      if( !paste("package",requiredPackage, sep=":") %in% search() ) attachNamespace(requiredPackage)
  }
  if(distribution=="gompertz"){
    requiredPackage<-"eha"
    if( !requireNamespace(requiredPackage, quietly = TRUE) ) 
      stop("Package \"",requiredPackage,"\" needed for option distribution=", distribution, ". Please install it.",
           call. = FALSE)
    else
      if( !paste("package",requiredPackage, sep=":") %in% search() ) attachNamespace(requiredPackage)
  }
  if(distribution=="pert"){
    requiredPackage<-"mc2d"
    if( !requireNamespace(requiredPackage, quietly = TRUE) ) 
      stop("Package \"",requiredPackage,"\" needed for option distribution=", distribution, ". Please install it.",
           call. = FALSE)
    else
      if( !paste("package",requiredPackage, sep=":") %in% search() ) attachNamespace(requiredPackage)
  }
  if(distribution=="tnorm"){
    requiredPackage<-"msm"
    if( !requireNamespace(requiredPackage, quietly = TRUE) ) 
      stop("Package \"",requiredPackage,"\" needed for option distribution=", distribution, ". Please install it.",
           call. = FALSE)
    else
      if( !paste("package",requiredPackage, sep=":") %in% search() ) attachNamespace(requiredPackage)
  }
  ## Consistency of arguments:
  if ( length(percentiles) != length(quantiles) )
    stop("length(percentiles) != length(quantiles)" )
  if( any( percentiles < 0 || percentiles > 1 ) )
    stop( "All elements of \"percentiles\" must lie between 0 and 1!")
  # Fit the distribution to the given quantiles:
  capture.output(dists<-try(rriskDistributions::rriskFitdist.perc(p=percentiles,
                                                                  q=quantiles,
                                                                  show.output=TRUE,
                                                                  tolConv=tolConv,
                                                                  fit.weights=fit.weights)),
                 file=if(verbosity>1) stdout() else NULL)
  if( inherits(dists, "try-error") )
    stop(dists$message)
  if(length(dists)==1 && is.na(dists))
    stop("No distribution could be fitted at absolute convergence tolerance tolConv=", tolConv, ".")
  ## Get the types of distributions that could be fitted:
  possible_dists<-colnames(dists$results[,3:ncol(dists$results)])[which(!is.na(dists$results[1,3:ncol(dists$results)]))]
  # Generate the random numbers according to the distribution type:
  if( match(distribution, possible_dists, nomatch=0) ){
    if(distribution=="norm") 
      x<-rnorm(n=n, 
               mean=dists$results[1,"norm"],
               sd=dists$results[2,"norm"])
    else if(distribution=="beta") 
      x<-rbeta(n=n,
               shape1=dists$results[1,"beta"],
               shape2=dists$results[2,"beta"])
    else if(distribution=="cauchy") 
      x<-rcauchy(n=n,
                 location=dists$results[1,"cauchy"],
                 scale=dists$results[2,"cauchy"])
    else if(distribution=="logis") 
      x<-rlogis(n=n,
                location=dists$results[1,"logis"],
                scale=dists$results[2,"logis"])
    else if(distribution=="t") 
      x<-rt(n=n,
            df=dists$results[1,"t"])
    else if(distribution=="chisq") 
      x<-rchisq(n=n,
                df=dists$results[1,"chisq"])
    else if(distribution=="chisqnc")
      x<-rchisq(n=n,
                df=dists$results[1,"chisqnc"],
                ncp=dists$results[2,"chisqnc"])
    else if(distribution=="exp") 
      x<-rexp(n=n,
              rate=dists$results[1,"exp"])
    else if(distribution=="f") 
      x<-rf(n=n,
            df1=dists$results[1,"f"],
            df2=dists$results[2,"f"])
    else if(distribution=="gamma") 
      x<-rgamma(n=n,
                shape=dists$results[1,"gamma"],
                rate=dists$results[2,"gamma"])
    else if(distribution=="lnorm") 
      x<-rlnorm(n=n,
                meanlog=dists$results[1,"lnorm"],
                sdlog=dists$results[2,"lnorm"])
    else if(distribution=="unif") {
      #unif needs exactly 2 quantiles (can't handle 3 or more)
      x<-runif(n=n,
               min=dists$results[1,"unif"],
               max=dists$results[2,"unif"])
    } else if(distribution=="weibull") 
      x<-rweibull(n=n,
                  shape=dists$results[1,"weibull"],
                  scale=dists$results[2,"weibull"])
    else if(distribution=="triang")
      x<-mc2d::rtriang(n=n,
                       min=dists$results[1,"triang"],
                       mode=dists$results[2,"triang"],
                       max=dists$results[3,"triang"])
    else if(distribution=="gompertz") 
      x<-eha::rgompertz(n=n,
                        shape=dists$results[1,"gompertz"],
                        scale=dists$results[2,"gompertz"])
    else if(distribution=="pert"){
      #pert needs 4 or more quantiles
      x<-mc2d::rpert(n=n,
                     min=dists$results[1,"pert"],
                     mode=dists$results[2,"pert"],
                     max=dists$results[3,"pert"],
                     shape=dists$results[4,"pert"])
    }
    else if(distribution=="tnorm") {
      #tnorm needs 4 or more quantiles
      x<-msm::rtnorm(n=n,
                     mean=dists$results[1,"tnorm"],
                     sd=dists$results[2,"tnorm"],
                     lower=dists$results[3,"tnorm"],
                     upper=dists$results[4,"tnorm"])
    } else
      stop("\"", distribution, "\" is not a valid distribution type.")
  } else {
    stop("\"", distribution, "\" distribution could not be fitted. One of the following should work:", 
         paste(possible_dists,collapse=", "))
  }
  # Check postcondition (i.e. goodness of fit):
  quantiles_calc <- quantile(x=x, probs=percentiles)
  percentiles_calc<-vapply(X=quantiles, 
                           FUN.VALUE=percentiles[[1]],
                           FUN=function(x_) length(x[x<=x_])/n
  )
  for( j in seq(along=quantiles) ){
        scale <- if( abs(quantiles[[j]]) > 0 ) quantiles[[j]] else NULL
        if( !isTRUE( msg<-all.equal(quantiles[[j]], quantiles_calc[[j]],  scale=scale, tolerance=relativeTolerance) ) ){
          warning("Fitted value of ", 100*percentiles[[j]], "%-quantile: ", quantiles_calc[[j]], "\n  ",
                  "Target value of ", 100*percentiles[[j]], "%-quantile: ", quantiles[[j]],   "\n  ",
                  "Fitted cumulative probability at value ", quantiles[[j]], " : ", percentiles_calc[[j]], "\n  ",
                  "Target cumulative probability at value ", quantiles[[j]], " : ", percentiles[[j]], "\n  ",
                  msg)
        }    
  }
  # Return sampled distribution if it could be achieved, NA otherwise:
  x
}
