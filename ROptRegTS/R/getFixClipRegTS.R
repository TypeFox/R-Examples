###############################################################################
## optimal clipping bound for finite-sample under-/overshoot risk
###############################################################################
setMethod("getFixClipRegTS", signature(clip = "numeric", 
                                       ErrorDistr = "Norm",
                                       Regressor = "UnivariateDistribution",
                                       risk = "fiUnOvShoot", 
                                       neighbor = "ContNeighborhood"),
    function(clip, ErrorDistr, Regressor, risk, neighbor){
        c0fct <- function(x, c0, tau){
            if(x == 0) return(0)

            return(exp(-2*c0*tau)*pnorm(tau*abs(x) - c0/abs(x)) 
                   - pnorm(-tau*abs(x)-c0/abs(x)))
        }
        rad <- neighbor@radius
        tau <- risk@width
        
        return(rad/(1-rad) - E(Regressor, c0fct, c0 = clip, tau = tau))
    })
setMethod("getFixClipRegTS", signature(clip = "numeric", 
                                       ErrorDistr = "Norm",
                                       Regressor = "UnivariateDistribution",
                                       risk = "fiUnOvShoot", 
                                       neighbor = "TotalVarNeighborhood"),
    function(clip, ErrorDistr, Regressor, risk, neighbor){
        c0fct <- function(x, c0, tau){
            if(x == 0) return(0)

            return(exp(-2*c0*tau)*pnorm(tau*abs(x) - c0/abs(x)) 
                   - pnorm(-tau*abs(x)-c0/abs(x)))
        }
        tau <- risk@width
        
        return((1 + exp(-2*clip*tau))*neighbor@radius 
               - E(Regressor, c0fct, c0 = clip, tau = tau))
    })
setMethod("getFixClipRegTS", signature(clip = "numeric", 
                                       ErrorDistr = "Norm",
                                       Regressor = "numeric",
                                       risk = "fiUnOvShoot", 
                                       neighbor = "CondContNeighborhood"),
    function(clip, ErrorDistr, Regressor, risk, neighbor){
        x <- Regressor
        rad <- neighbor@radiusCurve(x)
        tau <- risk@width
        res <- exp(-2*clip*abs(x)*tau)*pnorm(tau*abs(x) - clip) - pnorm(-tau*abs(x)-clip)
        
        return(rad/(1-rad) - res)
    })
setMethod("getFixClipRegTS", signature(clip = "numeric", 
                                       ErrorDistr = "Norm",
                                       Regressor = "numeric",
                                       risk = "fiUnOvShoot", 
                                       neighbor = "CondTotalVarNeighborhood"),
    function(clip, ErrorDistr, Regressor, risk, neighbor){
        x <- Regressor
        tau <- risk@width
        res <- (1 + exp(-2*tau*abs(x)*clip))*neighbor@radiusCurve(x)
        res <- res - exp(-2*clip*abs(x)*tau)*pnorm(tau*abs(x) - clip) + pnorm(-tau*abs(x)-clip)
        
        return(res)
    })
