inv_glob <- function (eta, type = "g", der = F) 
{
    D = NULL
    eta = as.vector(eta)
    if (length(eta) == 1) {
        ep = exp(eta)
        if(ep==Inf) p = c(0,1)
        else if(ep==-Inf) p = c(1,0)
        else p = c(1, ep)/(1 + ep)
        out = list(p = p, D = D)
    }
    else {
        if (length(eta) == 1) {
            p = exp(eta)
            p = p/(1 + p)
            p = c(1 - p, p)
            out = list(p = p, D = D)
        }
        else {
            if (type == "g") {
                eeta = exp(eta)
                pc = c(0, 1/(1 + eeta), 1)
                if(any(abs(eeta)==Inf)){
	                pc[eeta==Inf] = 1
	                pc[eeta==-Inf] = -1
                } 
                p = diff(pc)
                if (der) {
                  zz = rep(0, length(eta))
                  D = diff(rbind(zz, diag(-eeta/(1 + eeta)^2), 
                    zz))
                }
            }
            if (type == "l") {
                pc = exp(c(0, cumsum(eta)))
                p = pc/sum(pc)
            }
            out = list(p = p, D = D)
        }
    }
}