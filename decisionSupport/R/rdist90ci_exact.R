#
# file: rdist90ci_exact.R
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
# rdist90ci_exact(distribution, n, lower, upper)
##############################################################################################
#' 90\%-confidence interval based univariate random number generation (by exact parameter 
#' calculation).
#' 
#' This function generates random numbers for a set of univariate parametric distributions from  
#' given 90\% confidence interval.  Internally, this is achieved by exact, i.e. analytic, calculation
#' of the parameters for the individual distribution from the given 90\% confidence interval.
#' @param distribution \code{character}; A character string that defines the univariate distribution
#'  to be randomly sampled. For possible options cf. section Details.
#' @param n Number of generated observations.
#' @param lower \code{numeric}; lower bound of the 90\% confidence interval.
#' @param upper \code{numeric}; upper bound of the 90\% confidence interval.
#' @details
#'   The following table shows the available distributions and their identification 
#'  (option: \code{distribution}) as a character string:
#'  \tabular{lll}{
#'  \bold{\code{distribution}} \tab \bold{Distribution Name}      \tab \bold{Requirements}\cr
#'  \code{"const"}             \tab Deterministic case            \tab \code{lower == upper}\cr
#'  \code{"norm"}              \tab \link{Normal}                 \tab \code{lower < upper} \cr
#  \code{"pois"}              \tab ToDo \cr 
#  \code{"binom"}             \tab ToDo \cr 
#'  \code{"lnorm"}             \tab \link[=Lognormal]{Log Normal} \tab \code{0 < lower < upper} \cr
#  \code{"lnorm_lim2"}        \tab ToDo \cr
#'  \code{"unif"}              \tab \link{Uniform}                \tab \code{lower < upper}
#'  }
#'  \subsection{Parameter formulae}{
#'    We use the notation: \eqn{l}\code{=lower} and \eqn{u}=\code{upper}; 
#'    \eqn{\Phi} is the cumulative distribution function of the standard normal distribution and 
#'    \eqn{\Phi^{-1}}{\Phi^(-1)} its inverse, which is the quantile function of the standard normal 
#'    distribution.
#'    \describe{
#'      \item{\code{distribution="norm":}}{The formulae for \eqn{\mu} and \eqn{\sigma}, viz. the 
#'      mean and standard deviation, respectively, of the normal distribution are
#'      \eqn{\mu=\frac{l+u}{2}}{\mu=(l+u)/2} and
#'      \eqn{\sigma=\frac{\mu - l}{\Phi^{-1}(0.95)}}{\sigma=(\mu - l)/\Phi^(-1)(0.95)}.
#'      }
#'      \item{\code{distribution="unif":}}{For the minimum \eqn{a} and 
#'    maximum \eqn{b} of the uniform distribution \eqn{U_{[a,b]}}{U([a,b])} it holds that
#'                                          \eqn{a = l - 0.05 (u - l)
#'                                              }{a = l - 0.05 (u - l)} and 
#'                                          \eqn{b= u + 0.05 (u - l)
#'                                              }{b = u + 0.05 (u - l)}.
#'      }
#'       \item{\code{distribution="lnorm":}}{The density of the log normal distribution is 
#'             \eqn{ f(x) = \frac{1}{ \sqrt{2 \pi} \sigma x } \exp( - \frac{( \ln(x) - \mu )^2 %
#'             }{
#'             2 \sigma^2}) }{f(x)=1/((2 \pi)^(1/2)\sigma x) exp(-1/2(((ln(x)-\mu)/\sigma)^2))
#'             } for \eqn{x > 0} and \eqn{f(x) = 0} otherwise.
#'             Its parameters are determined by the confidence interval via 
#'             \eqn{\mu = \frac{\ln(l) + \ln(u)}{2}
#'             }{
#'             \mu = (ln(l)+ln(u))/2
#'             } and 
#'             \eqn{\sigma = \frac{1}{\Phi^{-1}(0.95)} ( \mu - \ln(l) )
#'             }{
#'             \sigma =  (\mu-ln(l))/\Phi^(-1)(0.95)
#'             }. Note the correspondence to the formula for the normal distribution.
#'       }
#'    }
#'  }
#' @return A numeric vector of length \code{n} with the sampled values according to the chosen 
#'   distribution.
#'   
#'   In case of \code{distribution="const"}, viz. the deterministic case, the function returns: 
#'   \code{rep(lower, n).}
#' @examples
#' # Generate uniformly distributed random numbers:
#' lower=3
#' upper=6
#' hist(r<-rdist90ci_exact(distribution="unif", n=10000, lower=lower, upper=upper),breaks=100)
#' print(quantile(x=r, probs=c(0.05,0.95)))
#' print(summary(r))
#' 
#' # Generate log normal distributed random numbers:
#' hist(r<-rdist90ci_exact(distribution="lnorm", n=10000, lower=lower, upper=upper),breaks=100)
#' print(quantile(x=r, probs=c(0.05,0.95)))
#' print(summary(r))
#' @export
rdist90ci_exact <- function(distribution, n, lower, upper){
  # Check preconditions
  if ( is.null(lower) || is.null(upper) || is.na(lower) || is.na(upper) )
    stop("lower and upper value of the 90%-confidence interval must be given.")
  # Prepare input variable: types
  lower<-as.numeric(lower)
  upper<-as.numeric(upper)
  
  # Create output vector for the random numbers to be generated
  x<-vector(length=n)
  # Generate the random numbers according to the distribution type:
  if(distribution=="const") {
    if( isTRUE(all.equal(lower,upper)) ) 
      x<-rep(lower,n)
    else
      stop("lower: ", lower, " is not equal to upper: ", upper)
  }
  else if (lower >= upper)
    stop("lower >= upper")
  else if(distribution=="norm") {    
    # 95%-critical value of standard normal distribution (c_0.95=1.645):
    c_0.95=qnorm(0.95)
    x<-rnorm(n=n,
             mean=mean(c(lower,upper)),
             sd=(mean(c(lower,upper))-lower)/c_0.95)
  }
  #   else if(distribution=="pois")    
  #     x<-rpois(n, mean(c(lower,upper))) 
  #   else if(distribution=="binom")   
  #     x<-rbinom(n,1,lower)  # ToDo: Why size=1? Why prob=lower?
  else if(distribution=="unif"){ 
    x<-runif(n=n, 
             min=lower-(upper-lower)*0.05, 
             max=upper+(upper-lower)*0.05)
    #old (wrong):    x<-runif(n,lower,upper)
  }
  else if(distribution=="lnorm"){ 
    if ( lower <= 0 )
      stop("lower <= 0")
    # 95%-critical value of standard normal distribution (c_0.95=1.645):
    c_0.95=qnorm(0.95)
    # Mean of defining normal distribution:
    meanlog<-mean( c(log(lower),log(upper)) )
    # Standard deviation of defining normal distribution:
    sdlog<-( meanlog - log(lower) )/c_0.95
    # Generate the random numbers:
    x<-rlnorm(n=n,
              meanlog = meanlog, 
              sdlog = sdlog)
    #     #old (wrong): 
    #x<-rlnorm(n,meanlog=mean(log(upper),log(lower)),sdlog=(mean(c(log(upper),log(lower)))-log(lower))/c_0.95) 
    # this is right (it was a typing error) and is the same as the new implementation:
    #x<-rlnorm(n,meanlog=mean(c(log(upper),log(lower))),sdlog=(mean(c(log(upper),log(lower)))-log(lower))/c_0.95) 
  }
  #   else if(distribution=="lnorm_lim2") {
  #     temp<-rlnorm(n,meanlog=mean(log(upper),log(lower)),sdlog=(mean(c(log(upper),log(lower)))-log(lower))/c_0.95)
  #     temp[which(temp>2*upper)]<-2*upper
  #     x<-temp
  #   }
  else
    stop("\"", distribution, "\" is not a valid distribution type.")
  #Return
  x
}
