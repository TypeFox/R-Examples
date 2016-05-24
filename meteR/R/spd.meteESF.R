
#' @title Species Power Distribution 
#'  
#' @description \code{spd.meteESF} calculates the distribution nu(e | N0, S0, E0), 
#' the distribution of average metabolic rates across for each species in a commmunity 
#' @details
#' See examples.
#' 
#' @param x an object of class meteESF. 
#' @keywords lagrange multiplier, METE, MaxEnt, ecosystem structure function
#' @export
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(spp=arth$spp,
#'                abund=arth$count,
#'                power=arth$mass^(.75),
#'                minE=min(arth$mass^(.75)))
#' spd1 <- spd(esf1)
#' 
#' @return An object of class \code{meteDist}. The object contains a list with the following elements.
#' \describe{
#'   \item{data}{The data used to construct the prediction}
#'   \item{d}{density funciton}
#'   \item{p}{cumulative density function}
#'   \item{q}{quantile funtion}
#'   \item{r}{random number generator}
#'   \item{La}{Vector of Lagrange multipliers}
#'   \item{state.var}{State variables used to constrain entropy maximization}
#'   \item{type}{Specifies the type of distribution is 'sad'}
#' }
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso meteDist, sad.meteESF, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.

spd <- function(x) {
  UseMethod('spd')
}

#' @rdname spd
#' @export 

spd.meteESF <- function(x) {
  if(is.na(x$state.var[3])) stop('must provide metabolic rate data or E0 to calculate power distributions')
  
  dat <- x$data$e
  
  if(is.null(dat)) {
    X <- NULL
  } else {
    X <- sort(tapply(dat, x$data$s, mean), decreasing=TRUE)
  }
  
  this.eq <- function(epsilon, log=FALSE) {
    out <- meteNu(epsilon, la1=x$La[1], 
                  la2=x$La[2], Z=x$Z,
                  S0=x$state.var[1], 
                  N0=x$state.var[2], 
                  E0=x$state.var[3])
    
    if(log) out <- log(out)
    
    return(out)
  }
  
  FUN <- distr::AbscontDistribution(d=this.eq,
                                    low1=1, low=1, up=x$state.var[3], up1=x$state.var[3])
  
  out <- list(type='spd', data=X, 
              d=this.eq, p=FUN@p, q=FUN@q, r=FUN@r,
              state.var=x$state.var, La=x$La)
  class(out) <- c('spd', 'meteDist')
  
  return(out)
}

##============================================================================
#' @title Equation of the PMF for the METE species metabolic rate distribution
#'
#' @description
#' \code{meteNu} is a low level function to calculate the value of
#' \eqn{\nu(e | N_0, S_0, E_0)} (the distribution of metabolic rates/power across all species in a commmunity) at the given value of \code{e}; vectorized in \code{e}.
#'
#' @details
#' Typically only used in \code{spd.meteESF} and not called by the user.
#' 
#' @param e the value (metabolic rate/power) at which to calculate \eqn{\Psi}
#' @param la1,la2 Lagrange multipliers
#' @param Z partition function
#' @param S0 Total number of species
#' @param N0 Total number of individuals
#' @param E0 Total metabolic rate
#' @keywords manip
#' @export
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(spp=arth$spp,
#'                 abund=arth$count,
#'                 power=arth$mass^(.75),
#'                 minE=min(arth$mass^(.75)))
#' meteNu(1:10,
#'        esf1$La[1],esf1$La[2],
#'        esf1$Z,esf1$state.var['S0'],
#'        esf1$state.var['N0'],
#'        esf1$state.var['E0'])
#' 
#' @return numeric vector of length equal to length of \code{e}
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#  @note other junk to mention
#' @seealso \code{spd.mete}
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
#  @aliases - a list of additional topic names that will be mapped to this documentation when the user looks them up from the command line.
#  @family - a family name. All functions that have the same family tag will be linked in the documentation.

meteNu <- function(e,la1,la2,Z,S0,N0,E0) {
  beta <- la1 + la2
  
  t1 <- 1/log(1/beta)
  t2 <- exp(-beta/(la2*(e-1)))
  t3 <- e-1
  
  out <- t1 * t2 / t3
  out[e==1] <- .Machine$double.eps

  return(out)
}

