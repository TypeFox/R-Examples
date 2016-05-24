###############################################################################
## optimal clipping bound for asymptotic MSE
###############################################################################
setMethod("getInfClipRegTS", signature(clip = "numeric", 
                                       ErrorL2deriv = "UnivariateDistribution",
                                       Regressor = "Distribution",
                                       risk = "asMSE", 
                                       neighbor = "Neighborhood"),
    function(clip, ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent){
        return(neighbor@radius^2*clip + getInfGammaRegTS(ErrorL2deriv = ErrorL2deriv, 
                                            Regressor = Regressor, risk = risk, 
                                            neighbor = neighbor, z.comp = z.comp, 
                                            stand = stand, cent = cent, clip = clip))
    })
setMethod("getInfClipRegTS", signature(clip = "numeric", 
                                       ErrorL2deriv = "UnivariateDistribution",
                                       Regressor = "Distribution",
                                       risk = "asMSE", 
                                       neighbor = "Av1CondTotalVarNeighborhood"),
    function(clip, ErrorL2deriv, Regressor, risk, neighbor, z.comp, stand, cent){
        if(!z.comp){
            cent <- function(x){-b/2}
            body(cent) <- substitute({-b/2}, list(b = clip))
            return(neighbor@radius^2*clip + 
                   getInfGammaRegTS(ErrorL2deriv = ErrorL2deriv, 
                                    Regressor = Regressor, risk = risk, 
                                    neighbor = neighbor, z.comp = z.comp, 
                                    stand = stand, cent = cent, clip = clip))
        }else{
            return(neighbor@radius^2*clip + 
                   getInfGammaRegTS(ErrorL2deriv = ErrorL2deriv, 
                                    Regressor = Regressor, risk = risk, 
                                    neighbor = neighbor, z.comp = z.comp, 
                                    stand = stand, cent = cent, clip = clip))       
        }
    })
setMethod("getInfClipRegTS", signature(clip = "numeric", 
                                       ErrorL2deriv = "EuclRandVariable",
                                       Regressor = "Distribution",
                                       risk = "asMSE", 
                                       neighbor = "Neighborhood"),
    function(clip, ErrorL2deriv, Regressor, risk, neighbor, ErrorDistr, stand, cent, trafo){
        return(neighbor@radius^2*clip + getInfGammaRegTS(ErrorL2deriv = ErrorL2deriv, 
                                            Regressor = Regressor, risk = risk, neighbor = neighbor, 
                                            ErrorDistr = ErrorDistr, stand = stand, cent = cent, 
                                            clip = clip))
    })
setMethod("getInfClipRegTS", signature(clip = "numeric", 
                                       ErrorL2deriv = "UnivariateDistribution",
                                       Regressor = "UnivariateDistribution",
                                       risk = "asUnOvShoot", 
                                       neighbor = "UncondNeighborhood"),
    function(clip, ErrorL2deriv, Regressor, risk, neighbor, z.comp, cent){
        if(!z.comp){
            return(neighbor@radius/risk@width + 
                   getInfGammaRegTS(ErrorL2deriv = ErrorL2deriv, 
                                    Regressor = Regressor, risk = risk, 
                                    neighbor = neighbor, cent = -clip/2, clip = clip))
        }else{
            return(neighbor@radius/risk@width + 
                   getInfGammaRegTS(ErrorL2deriv = ErrorL2deriv, 
                                    Regressor = Regressor, risk = risk, 
                                    neighbor = neighbor, cent = cent, clip = clip))       
        }
    })
setMethod("getInfClipRegTS", signature(clip = "numeric", 
                                       ErrorL2deriv = "UnivariateDistribution",
                                       Regressor = "numeric",
                                       risk = "asUnOvShoot", 
                                       neighbor = "CondNeighborhood"),
    function(clip, ErrorL2deriv, Regressor, risk, neighbor){
        return(neighbor@radiusCurve(Regressor)/risk@width + 
               getInfGammaRegTS(ErrorL2deriv = ErrorL2deriv, 
                                Regressor = Regressor, risk = risk, 
                                neighbor = neighbor, clip = clip))
    })
