###############################################################################
## optimal radius for given clipping bound for asymptotic MSE
###############################################################################

setMethod("getInfRad", signature(clip = "numeric",
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asMSE",
                                  neighbor = "ContNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype,
             cent, symm, trafo){
        gamm <- getInfGamma(L2deriv = L2deriv, risk = risk, neighbor = neighbor,
                            biastype = biastype, cent = cent, clip = clip)
        return((-gamm/clip)^.5)
    })

setMethod("getInfRad", signature(clip = "numeric",
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asMSE",
                                  neighbor = "TotalVarNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype,
             cent, symm, trafo){
        gamm <- getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk,
                            neighbor = neighbor, biastype = biastype,
                            cent = if(symm) -clip/2 else cent , clip = clip)
        return((-gamm/clip)^.5)
    })

setMethod("getInfRad", signature(clip = "numeric",
                                  L2deriv = "EuclRandVariable",
                                  risk = "asMSE",
                                  neighbor = "UncondNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype,
             Distr, stand, cent, trafo){
        gamm <- getInfGamma(L2deriv = L2deriv, risk = risk, neighbor = neighbor,
                            biastype = biastype, Distr = Distr, stand = stand,
                            cent = cent, clip = clip)
        return((-gamm/clip)^.5)
    })

###############################################################################
## optimal radius for given clipping bound for asymptotic L1risk
###############################################################################

setMethod("getInfRad", signature(clip = "numeric",
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asL1",
                                  neighbor = "ContNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype,
             cent, symm, trafo){
        s <- getInfV(L2deriv, neighbor, biastype, clip, cent, stand=1)
        gamm <- getInfGamma(L2deriv = L2deriv, risk = risk, neighbor = neighbor,
                            biastype = biastype, cent = cent, clip = clip)
        solvfct <- function(r){
            w <- r * clip / s^.5
            dp <- 2*dnorm(w)
            pp <- 2*pnorm(w)-1
            lhs <- s^.5*r*pp/dp
            return(lhs + gamm)
        }
        r <- try(uniroot(solvfct, lower=1e-5, upper = 10)$root,silent=TRUE)
        if(is(r, "try-error")) return(NA)
        return(r)
    })

setMethod("getInfRad", signature(clip = "numeric",
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asL1",
                                  neighbor = "TotalVarNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype,
             cent, symm, trafo){
        s <- getInfV(L2deriv, neighbor, biastype, clip, cent, stand=1)
        gamm <- getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk,
                               neighbor = neighbor, biastype = biastype,
                               cent = if(symm) -clip/2 else cent , clip = clip)
        solvfct <- function(r){
            w <- r * clip / s^.5
            dp <- 2*dnorm(w)
            pp <- 2*pnorm(w)-1
            lhs <- s^.5*r*pp/dp
            return(lhs + gamm)
        }
        r <- try(uniroot(solvfct, lower=1e-5, upper = 10)$root,silent=TRUE)
        if(is(r, "try-error")) return(NA)
        return(r)
    })


###############################################################################
## optimal radius for given clipping bound for asymptotic L4 risk
###############################################################################

setMethod("getInfRad", signature(clip = "numeric",
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asL4",
                                  neighbor = "ContNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype,
             cent, symm, trafo){
        s <- getInfV(L2deriv, neighbor, biastype, clip, cent, stand=1)
        gamm <- getInfGamma(L2deriv = L2deriv, risk = risk, neighbor = neighbor,
                            biastype = biastype, cent = cent, clip = clip)
        solvfct <- function(r){
           mse <- r^2 *clip^2 + s
           mse4 <- (r^2 *clip^2/3 + s)/mse
           r^2*clip*mse4 + gamm }
        r <- try(uniroot(solvfct, lower=1e-5, upper = 10)$root,silent=TRUE)
        if(is(r, "try-error")) return(NA)
        return(r)
    })

setMethod("getInfRad", signature(clip = "numeric",
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asL4",
                                  neighbor = "TotalVarNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype,
             cent, symm, trafo){
        s <- getInfV(L2deriv, neighbor, biastype, clip, cent, stand=1)
        gamm <- getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk,
                            neighbor = neighbor, biastype = biastype,
                            cent = if(symm) -clip/2 else cent , clip = clip)
        solvfct <- function(r){
           mse <- r^2 *clip^2 + s
           mse4 <- (r^2 *clip^2/3 + s)/mse
           r^2*clip*mse4 + gamm }
        r <- try(uniroot(solvfct, lower=1e-5, upper = 10)$root,silent=TRUE)
        if(is(r, "try-error")) return(NA)
        return(r)
    })



###############################################################################
## optimal radius for given clipping bound for asymptotic under-/overshoot risk
###############################################################################
setMethod("getInfRad", signature(clip = "numeric",
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asUnOvShoot",
                                  neighbor = "UncondNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype,
             cent, symm, trafo){
        gamm <- getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk,
                            neighbor = neighbor, biastype = biastype,
                            cent = if(symm) -clip/2 else cent , clip = clip)
        return( -gamm*risk@width)
    })

###############################################################################
## optimal radius for given clipping bound for asymptotic semivariance
###############################################################################
setMethod("getInfRad", signature(clip = "numeric",
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asSemivar",
                                  neighbor = "ContNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype, cent,  symm, trafo){
        biastype <- if(sign(risk)==1) positiveBias() else negativeBias()
        z0 <- getInfCent(L2deriv = L2deriv, risk = risk, neighbor = neighbor,
                         biastype = biastype,
                         clip = max(clip, 1e-4), cent = 0, trafo = trafo,
                         symm = symm, tol.z = 1e-6)

        ga <- getInfGamma(L2deriv = L2deriv, risk = risk, neighbor = neighbor,
                          biastype = biastype, cent = cent, clip = clip)


        solvfct <- function(r){
            si <- if (sign(risk)>0) 1 else -1
            v0 <- E(L2deriv, function(x) pmin( x-z0,  si*clip)^2 )
            s0 <- sqrt(v0)
            sv <- r * clip / s0
            er <- r^2 * clip + r * s0 * dnorm(sv) / pnorm(sv) + ga
            return(er)
        }
        r <- try(uniroot(solvfct, lower=1e-5, upper = 10)$root,silent=TRUE)
        if(is(r, "try-error")) return(NA)
        return(r)
     })
