###############################################################################
## optimal clipping bound for asymptotic MSE
###############################################################################

setMethod("getInfClip", signature(clip = "numeric", 
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asMSE", 
                                  neighbor = "ContNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype, 
             cent, symm, trafo){
        return(neighbor@radius^2*clip + 
               getInfGamma(L2deriv = L2deriv, risk = risk, 
                           neighbor = neighbor, biastype = biastype, cent = cent, clip = clip))
    })

setMethod("getInfClip", signature(clip = "numeric", 
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asMSE", 
                                  neighbor = "TotalVarNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype, 
             cent, symm, trafo){
        if(symm){
            return(neighbor@radius^2*clip + 
                   getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk, 
                               neighbor = neighbor, biastype = biastype, cent = -clip/2, clip = clip))
        }else{
            return(neighbor@radius^2*clip + 
                   getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk, 
                               neighbor = neighbor, biastype = biastype, cent = cent, clip = clip))
        }
    })

setMethod("getInfClip", signature(clip = "numeric", 
                                  L2deriv = "EuclRandVariable",
                                  risk = "asMSE", 
                                  neighbor = "UncondNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype, 
             Distr, stand, cent, trafo){
        return(neighbor@radius^2*clip +
                getInfGamma(L2deriv = L2deriv, risk = risk, neighbor = neighbor, 
                            biastype = biastype, Distr = Distr, stand = stand, 
                            cent = cent, clip = clip))
    })

###############################################################################
## optimal clipping bound for asymptotic L1risk
###############################################################################

setMethod("getInfClip", signature(clip = "numeric", 
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asL1", 
                                  neighbor = "ContNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype, 
             cent, symm, trafo){
        s <- getInfV(L2deriv, neighbor, biastype, clip, cent, stand=1)
        r <- neighbor@radius
        w <- r * clip / s^.5
        dp <- 2*dnorm(w)
        pp <- 2*pnorm(w)-1          
        return(s^.5*r*pp/dp + 
               getInfGamma(L2deriv = L2deriv, risk = risk, 
                           neighbor = neighbor, biastype = biastype, cent = cent, clip = clip))
    })

setMethod("getInfClip", signature(clip = "numeric", 
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asL1", 
                                  neighbor = "TotalVarNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype, 
             cent, symm, trafo){
        s <- getInfV(L2deriv, neighbor, biastype, clip, cent, stand=1)
        r <- neighbor@radius
        w <- r * clip / s^.5
        dp <- 2*dnorm(w)
        pp <- 2*pnorm(w)-1
        lhs <- s^.5*r*pp/dp
        if(symm){
            return(lhs + 
                   getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk, 
                               neighbor = neighbor, biastype = biastype, cent = -clip/2, clip = clip))
        }else{
            return(lhs + 
                   getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk, 
                               neighbor = neighbor, biastype = biastype, cent = cent, clip = clip))
        }
    })


###############################################################################
## optimal clipping bound for asymptotic L4 risk
###############################################################################

setMethod("getInfClip", signature(clip = "numeric", 
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asL4", 
                                  neighbor = "ContNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype, 
             cent, symm, trafo){
        s <- getInfV(L2deriv, neighbor, biastype, clip, cent, stand=1)
        r <- neighbor@radius
        mse <- r^2 *clip^2 + s
        mse4 <- (r^2 *clip^2/3 + s)/mse
        return(r^2*clip*mse4 + 
               getInfGamma(L2deriv = L2deriv, risk = risk, 
                           neighbor = neighbor, biastype = biastype, cent = cent, clip = clip))
    })

setMethod("getInfClip", signature(clip = "numeric", 
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asL4", 
                                  neighbor = "TotalVarNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype, 
             cent, symm, trafo){
        s <- getInfV(L2deriv, neighbor, biastype, clip, cent, stand=1)
        r <- neighbor@radius
        mse <- r^2 *clip^2 + s
        mse4 <- (r^2 *clip^2/3 + s)/mse
        if(symm){
            return(r^2*clip*mse4 + 
                   getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk, 
                               neighbor = neighbor, biastype = biastype, cent = -clip/2, clip = clip))
        }else{
            return(r^2*clip*mse4 + 
                   getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk, 
                               neighbor = neighbor, biastype = biastype, cent = cent, clip = clip))
        }
    })



###############################################################################
## optimal clipping bound for asymptotic under-/overshoot risk
###############################################################################
setMethod("getInfClip", signature(clip = "numeric", 
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asUnOvShoot", 
                                  neighbor = "UncondNeighborhood"),
    function(clip, L2deriv, risk, neighbor, biastype, 
             cent, symm, trafo){
        if(symm){
            return(neighbor@radius/risk@width + 
                   getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk, 
                               neighbor = neighbor, biastype = biastype, cent = -clip/2, clip = clip))
        }else{
            return(neighbor@radius/risk@width + 
                   getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk, 
                               neighbor = neighbor, biastype = biastype, cent = cent, clip = clip))
        }
    })

###############################################################################
## optimal clipping bound for asymptotic semivariance
###############################################################################
setMethod("getInfClip", signature(clip = "numeric", 
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

        r <- neighbor@radius

        if (sign(risk)>0)
            v0 <- E(L2deriv, function(x) pmin( x-z0,  clip)^2 )
        else
            v0 <- E(L2deriv, function(x) pmax( x-z0, -clip)^2 )        

        s0 <- sqrt(v0)
        sv <- r * clip / s0

        er <- r^2 * clip + r * s0 * dnorm(sv) / pnorm(sv) + ga 
        return(er)
     })
