#' @title meteESF
#'  
#' @description \code{meteESF} Calculates the ``ecosystem structure
#' function'' \eqn{R(n,\epsilon)} which forms the core of the Maximum Entropy Theory of
#' Ecology
#'
#' @details
#' Uses either data or state variables to calculate the Ecosystem Structure 
#' Function (ESF). \code{power} nor \code{E0} need not be specified; if missing an arbitrarily
#' large value is assigned to E0 (N0*1e5) such that it will minimally affect 
#' estimation of Lagrange multipliers. Consider using sensitivity analysis to 
#' confirm this assumption. Examples show different ways of combining data and state
#' variables to specify constraints
#' 
#' @param spp A vector of species names 
#' @param abund A vector of abundances 
#' @param power A vector of metabolic rates 
#' @param S0 Total number of species
#' @param N0 Total number of individuals
#' @param E0 Total metabolic rate; defaults to N0*1e6 if not specified or 
#'                  calculated from \code{power} to allow one to fit models that 
#'                  do not depend on metabolic rates
#' @param minE Minimum possible metabolic rate
#' @keywords lagrange multiplier, METE, MaxEnt, ecosystem structure function
#' @export
#' 
#' @examples
#' ## case where complete data availible
#' esf1 <- meteESF(spp=arth$spp,
#'                 abund=arth$count,
#'                 power=arth$mass^(.75),
#'                 minE=min(arth$mass^(.75)))
#' esf1
#' 
#' ## excluding metabolic rate data
#' esf2 <- meteESF(spp=arth$spp,
#'                 abund=arth$count)
#' esf2
#' 
#' ## using state variables only
#' esf3 <- meteESF(S0=50, N0=500, E0=5000)
#' esf3
#' esf4 <- meteESF(S0=50, N0=500)
#' esf4
#' 
#' @return An object of class \code{meteESF} with elements
#' \describe{
#'   \item{\code{data}}{The data used to construct the ESF}
#'   \item{\code{emin}}{The minimum metabolic rate used to rescale metabolic rates}
#'   \item{\code{La}}{Vector of Lagrange multipliers}
#'   \item{\code{La.info}}{Termination information from optimization procedure}
#'   \item{\code{state.var}}{State variables used to constrain entropy maximization}
#'   \item{\code{Z}}{Normalization constant for ESF}
#' }
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso metePi
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.

meteESF <- function(spp, abund, power,
                    S0=NULL, N0=NULL, E0=NULL,
                    minE) {
    ## case where spp (and abund) provided
    if(!missing(spp) & !missing(abund)) {
        ## factors will mess things up
        spp <- as.character(spp)
        
        ## account for possible 0 abundances
        abund0 <- abund == 0
        spp <- spp[!abund0]
        abund <- abund[!abund0]
        if(!missing(power)) power <- power[!abund0]
        
        ## if power missing, set large for numeric approx
        if(missing(power)) {
            power <- rep(10, length(spp))
            power[1] <- 1 # do this so re-scaling by min still allows E0 big
            e.given <- FALSE # to determine if power should be returned as data
        } else {
            e.given <- TRUE
        }
        
        ## account for possible aggregation of individuals
        if(any(abund > 1)) {
            spp <- rep(spp, abund)
            power <- rep(power, abund)
            abund <- rep(1, length(spp))
        }
        
        ## if no theoretical minimum metabolic rate set, set to min of
        ## data
        if(missing(minE)) minE <- min(power)
        power <- power/minE
        
        ## combine data, potentially e too if given
        dats <- data.frame(s=spp, n=abund)
        if(e.given) dats$e <- power
        
        ## calculate state variables from data
        S0 <- length(unique(spp))
        N0 <- sum(abund)
        E0 <- sum(power*abund)
        
    } else if(is.null(E0)) {
        E0 <- N0*10^3 # make very large so it has no effect
        e.given <- FALSE
        dats <- NULL
    } else {
        e.given <- TRUE
        dats <- NULL
    }
    
    ## if E0 given but no minE set it here to 1
    if(missing(minE)) minE <- 1
    
    ## calculate ecosystem structure funciton
    thisESF <- .makeESF(s0=S0, n0=N0, e0=E0)
    thisESF$emin <- minE
    
    ## include data if applicable
    if(exists('dats')) thisESF$data <- dats
    
    ## to be output
    out <- thisESF
    
    ## remove E0 state variable if it was in fact not given
    if(!e.given) out$state.var['E0'] <- NA
    
    ## make returned object of class METE
    class(out) <- 'meteESF'
    
    return(out)
}

##=============================================================================
## helper fun to calculate computational costly parts of ESF: the lagrange multipliers and Z

.makeESF <- function(s0,n0,e0) {
    esf.par <- .mete.lambda(s0,n0,e0)
    esf.par.info <- esf.par[-1]
    esf.par <- esf.par[[1]]
    
    names(esf.par) <- c("la1","la2")
    
    thisZ <- .meteZ(esf.par["la1"],esf.par["la2"],s0,n0,e0)
    
    return(list(La=esf.par,La.info=esf.par.info,Z=thisZ,state.var=c(S0=s0,N0=n0,E0=e0)))
}


##===========================================================================
## helper fun to calculate lagrange multipliers
#' @importFrom stats nlm

.mete.lambda <- function(S0, N0, E0) {
    ## reasonable starting values
    init.la2 <- S0/(E0-N0)
    beta.guess <- 0.01
    
    init.beta <- nlm(function(b) {
        ## options set to surpress warning messages just for optimization
        orig.warn <- options(warn=-1)
        
        out <- (b*log(1/b) - 1/2^8)^2
        
        return(out)
    }, p=0.001)
    
    if(init.beta$code < 4) {	# was there some level of convergence?
        init.beta <- init.beta$estimate
    } else {
        init.beta <- beta.guess
    }
    
    init.la1 <- init.beta - init.la2
    
    ## the solution
    la.sol <- nleqslv::nleqslv(x=c(la1 = init.la1, la2 = init.la2), 
                               fn = .la.syst2, S0=S0,N0=N0,E0=E0)
    
    return(list(lambda=c(la1=la.sol$x[1], la2=la.sol$x[2]), 
                syst.vals=la.sol$fvec, converg=la.sol$termcd,
                mesage=la.sol$message, nFn.calc=la.sol$nfcnt, nJac.calc=la.sol$njcnt))
}


##==========================================================================
## system of equations needed to be solved for getting lagrange multipliers

.la.syst2 <- function(La, S0, N0, E0) {
    ## options set to surpress warning messages just for optimization
    orig.warn <- options(warn=-1)
    
    ## params
    b <- La[1] + La[2]
    s <- La[1] + E0*La[2]
    
    n <- 1:N0
    
    ## expressions
    g.bn <- exp(-b*n)
    g.sn <- exp(-s*n)
    
    univ.denom <- sum((g.bn - g.sn)/n)
    rhs.7.19.num <- sum(g.bn - g.sn)
    rhs.7.20.num <- sum(g.bn - E0*g.sn)
    
    ##  the two functions to solve
    f <- rep(NA, 2)
    f[1] <- rhs.7.19.num/univ.denom - N0/S0
    f[2] <- (1/La[2]) + rhs.7.20.num/univ.denom - E0/S0
    
    ## return warning behavior to original
    options(warn=orig.warn$warn)
    
    return(f)
}


##===========================================================================
## function to return Z, the normalizing constant of R as well as
## the simplifying parameters beta and sigma

.meteZ <- function(la1,la2,S0,N0,E0) {
    beta <- la1 + la2
    sigma <- la1 + E0*la2
    
    t1 <- S0/(la2*N0)
    t2 <- (exp(-beta) - exp(-beta*(N0+1)))/(1-exp(-beta))
    t3 <- (exp(-sigma) - exp(-sigma*(N0+1)))/(1-exp(-sigma))
    
    Z <- t1*(t2 - t3)
    
    return(Z)
}

##===========================================================================
#' @title predictESF
#'  
#' @description \code{predict} predicts the probabilities for given combinations of abundance and energ from the  ``ecosystem structure
#' function'' \eqn{R(n,\epsilon)} 
#'
#' @details
#' Uses a fitted object of class \code{meteESF} and user supplied values of abundance and power to predict values of the ESF
#' 
#' @param esf A fitted object of class \code{meteESF}
#' @param abund A vector of abundances 
#' @param power A vector of metabolic rates 
#' 
#' @keywords lagrange multiplier, METE, MaxEnt, ecosystem structure function
#' @export
#' 
#' @examples
#' ## case where complete data availible
#' esf1 <- meteESF(spp=arth$spp,
#'                 abund=arth$count,
#'                 power=arth$mass^(.75),
#'                 minE=min(arth$mass^(.75)))
#' predictESF(esf1,
#'            abund=c(10,3),
#'            power=c(.01,3))

#' 
#' @return a data.frame with abundance, power, and the predicted value of the ESF
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso meteESF
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.

predictESF=function(esf,abund,power){
  p=(1/esf$Z) * 
    exp(-1*esf$La[1]*abund) * 
    exp(-1*esf$La[1]*power)
  return(data.frame(abund=abund,
                    power=power,
                    ESF.prob=p))
}

