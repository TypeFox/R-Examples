
#' @title Individual Power Distribution 
#'  
#' @description \code{ipd.meteESF} calculates the distribution Psi(e | N0, S0, E0), 
#' the distribution of metabolic rates across all individuals in a commmunity 
#' @details
#' See examples.
#' 
#' @param x an object of class meteESF. 
#' @param ... additiona arguments to be passed to methods
#' 
#' @keywords lagrange multiplier, METE, MaxEnt, ecosystem structure function
#' @export
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(spp=arth$spp,
#'                abund=arth$count,
#'                power=arth$mass^(.75),
#'                minE=min(arth$mass^(.75)))
#' ipd1 <- ipd(esf1)
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

ipd <- function(x, ...) {
	UseMethod('ipd')
}

#' @rdname ipd
# @method ipd meteESF
# @S3method ipd meteESF

#' @export

#' @importFrom stats integrate
ipd.meteESF <- function(x,...) {
    if(is.na(x$state.var[3])) stop('must provide metabolic rate data or E0 to calculate power distributions')
    
    dat <- x$data$e
    if(is.null(dat)) {
        X <- NULL
    } else {
        X <- sort(dat, decreasing=TRUE)
    }
    
    this.eq <- function(epsilon, log=FALSE) {
        out <- metePsi(epsilon, la1=x$La[1], 
                       la2=x$La[2], Z=x$Z,
                       S0=x$state.var[1], 
                       N0=x$state.var[2], 
                       E0=x$state.var[3])
        
        if(log) out <- log(out)
        
        return(out)
    }

    this.p.eq <- Vectorize(function(epsilon, lower.tail=TRUE, log.p=FALSE) {
        out <- integrate(this.eq,lower=1, upper=epsilon)$value
        
        if(!lower.tail) out <- 1 - out
        
        if(log.p) out <- log(out)
        
        return(out)
    }, vectorize.args="epsilon")
    
    FUN <- distr::AbscontDistribution(d=this.eq, p=this.p.eq,
                                      low1=1, low=1, up=x$state.var[3], up1=x$state.var[3])
    
    out <- list(type='ipd', data=X, 
                d=this.eq, p=FUN@p, q=FUN@q, r=FUN@r,
                state.var=x$state.var, La=x$La)
    class(out) <- c('ipd', 'meteDist')
    
    return(out)
}

##============================================================================
#' @title Equation of the PMF for the METE individual metabolic rate distribution
#'
#' @description
#' \code{metePsi} is a low level function to calculate the value of
#' \eqn{\Psi(e | N_0, S_0, E_0)} (the distribution of metabolic rates/power across all individuals in a commmunity) at the given value of \code{e}; vectorized in \code{e}.
#'
#' @details
#' Typically only used in \code{ipd.meteESF} and not called by the user.
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
#' metePsi(1:10,
#'         esf1$La[1],esf1$La[2],
#'         esf1$Z,esf1$state.var['S0'],
#'         esf1$state.var['N0'],
#'         esf1$state.var['E0'])
#' 
#' @return numeric vector of length equal to length of \code{e}
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#  @note other junk to mention
#' @seealso \code{ipd.mete}
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
#  @aliases - a list of additional topic names that will be mapped to this documentation when the user looks them up from the command line.
#  @family - a family name. All functions that have the same family tag will be linked in the documentation.

metePsi <- function(e,la1,la2,Z,S0,N0,E0) {
    gamma <- la1 + e*la2
    
    t1 <- S0/(N0*Z)
    t2 <- exp(-gamma)/((1-exp(-gamma))^2)
    t3 <- exp(-gamma*N0)/(1-exp(-gamma))
    t4 <- N0 + exp(-gamma)/(1-exp(-gamma))
    
    return(t1*(t2 - t3*t4))
}

