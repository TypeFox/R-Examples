###############################################################################
## optimal clipping bound for finite-sample under-/overshoot risk
###############################################################################
setMethod("getFixClip", signature(clip = "numeric", 
                                  Distr = "Norm",
                                  risk = "fiUnOvShoot", 
                                  neighbor = "ContNeighborhood"),
    function(clip, Distr, risk, neighbor){
        rad <- neighbor@radius
        tau <- risk@width
        return(rad/(1-rad) + pnorm(-tau-clip) - exp(-2*clip*tau)*pnorm(tau-clip))
    })

setMethod("getFixClip", signature(clip = "numeric", 
                                  Distr = "Norm",
                                  risk = "fiUnOvShoot", 
                                  neighbor = "TotalVarNeighborhood"),
    function(clip, Distr, risk, neighbor){
        tau <- risk@width
        return((1 + exp(-2*clip*tau))*neighbor@radius + pnorm(-tau-clip) 
               - exp(-2*clip*tau)*pnorm(tau-clip))
    })
