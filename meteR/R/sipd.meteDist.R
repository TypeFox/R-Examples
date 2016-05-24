
#' @title Generic method to obtain the species-level individual power distribution (SIPD)
#'
#' @description Extract species level individual power distribution 
#' from ESF object and return object inheriting from meteDist. This distribution (Theta) describes the distribution of metabolic rates across the individuals of a species with n individuls
#'
#' @details 
#' If \code{n} is provided then only the theoretical prediction is returned (because 
#' data from multiple species could map to the same n). Thus if data and prediction are
#' desired use \code{sppID}.
#' 
#' 
#' @param x An object of class meteESF (i.e. the fitted distribution \eqn{R(n,e)})
#' @param sppID the name or index of the species of interest as listed in the \code{spp} argument passed to \code{meteESF}
#' @param n integer. Alternatively can extract METE prediction by indicating number of individuals in the species
#' @param ... arguments to be passed to methods
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
#' @rdname sipd
#' @export
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(spp=arth$spp,
#'                 abund=arth$count,
#'                 power=arth$mass^(.75),
#'                 minE=min(arth$mass^(.75)))
#' sipd1 <- sipd(esf1, sppID=5)
#' sipd1
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso sad.meteESF, ipd.meteESF, metePsi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
#' @family Theta

sipd <- function(x, ...) {
	UseMethod('sipd')
}

#' @rdname sipd
# @method sipd meteESF
# @S3method sipd meteESF
#' @export 

sipd.meteESF <- function(x, sppID, n,...) {
    if(is.na(x$state.var[3])) stop('must provide metabolic rate data or E0 to calculate power distributions')
    
    if(!missing(n)) { # n provided, only return theoretical (cause how deal with multiple spp with same n?)
        if(!missing(sppID)) warning('`sppID` ignored if `n` provided; to explore individual species use `sppID` instead')
        X <- NULL # data null if n provided
    } else if(!is.null(x$data)) { # data provided
        if(!missing(sppID)) {
            spp <- as.factor(x$data$s)
            if(is.numeric(sppID)) sppID <- levels(spp)[sppID]
            
            n <- sum(x$data$n[x$data$s == sppID])
            dat <- x$data$e[x$data$s == sppID]
            if(is.null(dat)) { # no energy data
                X <- NULL
            } else {
                X <- sort(dat, decreasing=TRUE)
            }
        } else {
            
            stop('must provide either `sppID` or `n`')
        }
    } else {
        stop('must provide `n` if `esf` has no empirical data')
    }
    
    this.eq <- function(epsilon, log=FALSE) {
        out <- meteTheta(epsilon, n=n, la2=x$La[2])
        if(log) out <- log(out)
        return(out)
    }
    
    this.p.eq <- function(epsilon, log=FALSE, lower.tail=TRUE) {
        out <- .meteThetaCum(epsilon, n, x$La[2], log)
        if(!lower.tail) out <- 1 - out
        
        return(out)
    }
    
    FUN <- distr::AbscontDistribution(d=this.eq, p=this.p.eq,
                                      low1=1, low=1, up=x$state.var[3], up1=x$state.var[3])
    
    out <- list(type='sipd', data=X,
                d=this.eq, p=this.p.eq, q=FUN@q, r=FUN@r,
                state.var=x$state.var, La=x$La)
    out$state.var['n'] <- n
    class(out) <- c('ipd', 'meteDist')
    
    return(out)
}

#================================================================
#' @title Equation of the PMF for the METE Intra-specific metabolic rate distribution
#'
#' @description Distribution of metabolic rates over individuals within a species of abundance n0 
#'
# @details  
# \deqn{
#    \Theta( \epsilon \mid n, S_{0}, N_{0}, E_{0} ) \approx \lambda_{2} n e^{- \lambda_[2] n (\epsilon -1)}
# }
#' 
#' @param e Metabolic rate
#' @param n Number of individuals in species
#' @param la2 Lagrange multiplier (lambda_2) as obtained from \code{meteESF}
#' 
#' @export
#' 
#' @return numeric vector of length equal to lengthof \code{e}
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso metePsi, ipd
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
#' @family Theta

meteTheta <- function(e, n, la2) {
    la2 * n * exp(-la2 * n * (e-1))
}

## integral of meteTheta
.meteThetaCum <- function(e, n, la2, log=FALSE) {
    out <- 1 - exp(-la2 * n * (e - 1))
    if(log) out <- log(out)
    
    return(out)
}
