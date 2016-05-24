###############################################################################
## optimal clipping bound for asymptotic MSE
###############################################################################
setMethod("getInfClip", signature(clip = "numeric", 
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asMSE", 
                                  neighbor = "ContNeighborhood"),
    function(clip, L2deriv, risk, neighbor, cent, symm, trafo){
        return(neighbor@radius^2*clip + 
               getInfGamma(L2deriv = L2deriv, risk = risk, 
                           neighbor = neighbor, cent = cent, clip = clip))
    })
setMethod("getInfClip", signature(clip = "numeric", 
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asMSE", 
                                  neighbor = "TotalVarNeighborhood"),
    function(clip, L2deriv, risk, neighbor, cent, symm, trafo){
        if(symm){
            return(neighbor@radius^2*clip + 
                   getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk, 
                               neighbor = neighbor, cent = -clip/2, clip = clip))
        }else{
            return(neighbor@radius^2*clip + 
                   getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk, 
                               neighbor = neighbor, cent = cent, clip = clip))
        }
    })
setMethod("getInfClip", signature(clip = "numeric", 
                                  L2deriv = "EuclRandVariable",
                                  risk = "asMSE", 
                                  neighbor = "ContNeighborhood"),
    function(clip, L2deriv, risk, neighbor, Distr, stand, cent, trafo){
        return(neighbor@radius^2*clip + 
                getInfGamma(L2deriv = L2deriv, risk = risk, neighbor = neighbor, 
                            Distr = Distr, stand = stand, cent = cent, clip = clip))
    })

###############################################################################
## optimal clipping bound for asymptotic under-/overshoot risk
###############################################################################
setMethod("getInfClip", signature(clip = "numeric", 
                                  L2deriv = "UnivariateDistribution",
                                  risk = "asUnOvShoot", 
                                  neighbor = "UncondNeighborhood"),
    function(clip, L2deriv, risk, neighbor, cent, symm, trafo){
        if(symm){
            return(neighbor@radius/risk@width + 
                   getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk, 
                               neighbor = neighbor, cent = -clip/2, clip = clip))
        }else{
            return(neighbor@radius/risk@width + 
                   getInfGamma(L2deriv = sign(as.vector(trafo))*L2deriv, risk = risk, 
                               neighbor = neighbor, cent = cent, clip = clip))
        }
    })
