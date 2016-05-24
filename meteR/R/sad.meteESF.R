
#' @title METE species abundance distribution
#'
#' @description
#' \code{sad.mete} returns the species abundance distribution
#' predicted by METE (\eqn{\Phi(n)})
#'
#' @details
#' See Examples.
#' 
#' @param x an object of class mete. 
#' @keywords lagrange multiplier, METE, MaxEnt, ecosystem structure 
#' function
#' @export
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(spp=arth$spp,
#'                 abund=arth$count,
#'                 power=arth$mass^(.75),
#'                 minE=min(arth$mass^(.75)))
#' sad1 <- sad(esf1)
#' sad1
#' sad1$r(20)
#' sad1$q(seq(0,1,length=10))
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
# @note
#' @seealso metePhi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

sad <- function(x) {
	UseMethod('sad')
}

#' @rdname sad
# @method sad meteESF
# @S3method sad meteESF
#' @export 

sad.meteESF <- function(x) {
    dat <- x$data$n
    otu <- x$data$s
    
    if(is.null(dat)) {
      X <- NULL
    } else {
      X <- sort(tapply(dat,otu,sum),decreasing=TRUE)
    }
    
    
    this.eq <- function(n, log=FALSE) {
      out <- metePhi(n=n,la1=x$La[1],la2=x$La[2],Z=x$Z,
                     S0=x$state.var[1],N0=x$state.var[2],
                     E0=ifelse(is.na(x$state.var[3]), x$state.var[2]*10^2, x$state.var[3]))
#      out=out/sum(out) # to address error below:  Error in distr::DiscreteDistribution(supp = 1:x$state.var["N0"], prob = this.eq(1:x$state.var["N0"])) : sum of 'prob' has to be (approximately) 1
      if(log) out <- log(out)
      return(out)
    }
    
    FUN <- distr::DiscreteDistribution(supp=1:x$state.var["N0"], 
                                prob=this.eq(1:x$state.var["N0"]))
    
    # rankFun <- qphi2rank(FUN@q,x$state.var["S0"])
    
    out <- list(type='sad', data=X, 
                d=this.eq, p=FUN@p, q=FUN@q, r=FUN@r,
                state.var=x$state.var, La=x$La)
    class(out) <- c('sad', 'meteDist')
    
    return(out)
}



#==============================================================================
#' @title Equation of the METE species abundance distribution
#'
#' @description
#' \code{metePhi} returns the species abundance distribution 
#' (Phi(n)) predicted by METE; vectorized in n
#'
#' @details
#' See Examples
#' 
#' @param n the value (number of individuals) at which to calculate \deqn{\Phi}
#' @param la1,la2 Lagrange multipliers
#' @param Z partition function
#' @param S0 Total number of species
#' @param N0 Total number of individuals
#' @param E0 Total metabolic rate
#' 
#' @keywords manip
#' @export
#' 
#' @examples
#' esf1=meteESF(spp=arth$spp,
#'               abund=arth$count,
#'               power=arth$mass^(.75),
#'               minE=min(arth$mass^(.75)))
#' metePhi(min(arth$mass^(.75)),
#'        esf1$La[1],esf1$La[2],
#'        esf1$Z,esf1$state.var['S0'],
#'        esf1$state.var['N0'],
#'        esf1$state.var['E0'])
#' 
#' @return numeric
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#  @note other junk to mention
#' @seealso \code{sad.mete}
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
#  @aliases - a list of additional topic names that will be mapped to this documentation when the user looks them up from the command line.
#  @family - a family name. All functions that have the same family tag will be linked in the documentation.


metePhi <- function(n, la1, la2, Z, S0, N0, E0) {
    if(missing(Z)) Z <- .meteZ(la1, la2, S0, N0, E0)
    
    beta <- la1 + la2
    sigma <- la1 + E0*la2
    
    return((exp(-beta*n) - exp(-sigma*n))/(la2*Z*n))
}
