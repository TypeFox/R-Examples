###############################################################################
## lower case radius
###############################################################################
setMethod("lowerCaseRadius", signature(L2Fam = "L2ParamFamily",
                                       neighbor = "ContNeighborhood",
                                       risk = "asMSE"),
    function(L2Fam, neighbor, risk){
        if(length(L2Fam@param) != 1) stop("not yet implemented")

        D1 <- L2Fam@distribution
        if(!is(D1, "DiscreteDistribution")) stop("not yet implemented")

        w0 <- options("warn")
        options(warn = -1)
        L2deriv <- L2Fam@L2derivDistr[[1]]        
        m <- q(L2deriv)(0.5)
        wsm <- d(L2deriv)(m)
        
        supp <- support(L2deriv)
        gg1 <- min(supp[supp > m] - m)
        gg2 <- max(supp[supp < m] - m)
        Int <- E(L2deriv, function(x, m){abs(x-m)}, m = m)

        if(wsm > 0){
            gg <- min(abs(supp[supp != m] - m))
            if(gg > 0){
                bet <- (2*p(L2deriv)(m) - wsm - 1)/wsm
                if((p(L2deriv)(m)-wsm == 0)|(1-p(L2deriv)(m) == 0)){
                    M <- (1 + abs(bet))/gg
                }else{
                    M <- max((1-bet)/gg1, (1+bet)/(-gg2))
                }
                rad <- sqrt(M*Int - (1-wsm) - bet^2*wsm)
                names(rad) <- "lower case radius"
                
                options(w0)
                return(rad)
            }else{
                options(w0)
                rad <- Inf
                names(rad) <- "lower case radius"
                return(rad)
            }
        }else{
            M <- 2/(gg1-gg2)
            rad <- sqrt(M*Int - 1)
            names(rad) <- "lower case radius"

            options(w0)
            return(rad)            
        }
    })
setMethod("lowerCaseRadius", signature(L2Fam = "L2ParamFamily",
                                       neighbor = "TotalVarNeighborhood",
                                       risk = "asMSE"),
    function(L2Fam, neighbor, risk){
        if(length(L2Fam@param) != 1) stop("not yet implemented")

        D1 <- L2Fam@distribution
        if(!is(D1, "DiscreteDistribution")) stop("not yet implemented")

        L2deriv <- L2Fam@L2derivDistr[[1]]        
        w0 <- options("warn")
        options(warn = -1)
        supp <- support(L2deriv)
        gg <- min(abs(supp[supp != 0]))
        if(gg > 0){
            gg1 <- min(supp[supp > 0])
            gg2 <- max(supp[supp < 0])
            ws0 <- 1 - d(L2deriv)(0)
            ws1 <- 1 - p(L2deriv)(0)
            ws2 <- p(L2deriv)(0) - d(L2deriv)(0)
            M <- 1/ws0*max(ws2/gg1, ws1/(-gg2))
            rad <- sqrt(M*(-m1df(L2deriv, 0)) - ws1*ws2/ws0)
            names(rad) <- "lower case radius"

            options(w0)
            return(rad)
        }else{
            rad <- Inf
            names(rad) <- "lower case radius"

            options(w0)
            return(rad)
        }
    })
