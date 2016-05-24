###############################################################################
## lower case radius
###############################################################################
setMethod("lowerCaseRadius", signature(L2Fam = "L2ParamFamily",
                                       neighbor = "ContNeighborhood",
                                       risk = "asMSE",
                                       biastype = "ANY"),
    function(L2Fam, neighbor, risk, biastype){
        if(length(L2Fam@param) != 1) stop("not yet implemented")

        D1 <- L2Fam@distribution
        if(!is(D1, "DiscreteDistribution")) stop("not yet implemented")

        w0 <- options("warn")
        on.exit(options(w0))
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
                
                return(rad)
            }else{
                rad <- Inf
                names(rad) <- "lower case radius"
                return(rad)
            }
        }else{
            M <- 2/(gg1-gg2)
            rad <- sqrt(M*Int - 1)
            names(rad) <- "lower case radius"

            return(rad)            
        }
    })
setMethod("lowerCaseRadius", signature(L2Fam = "L2ParamFamily",
                                       neighbor = "TotalVarNeighborhood",
                                       risk = "asMSE",
                                       biastype = "ANY"),
    function(L2Fam, neighbor, risk, biastype){
        if(length(L2Fam@param) != 1) stop("not yet implemented")

        D1 <- L2Fam@distribution
        if(!is(D1, "DiscreteDistribution")) stop("not yet implemented")

        L2deriv <- L2Fam@L2derivDistr[[1]]        
        w0 <- options("warn")
        on.exit(options(w0))
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
###############################################################################
# onesided and asymmetric terms
###############################################################################
setMethod("lowerCaseRadius", signature(L2Fam = "L2ParamFamily",
                                       neighbor = "ContNeighborhood",
                                       risk = "asMSE",
                                       biastype = "onesidedBias"),
    function(L2Fam, neighbor, risk, biastype){
        if(length(L2Fam@param) != 1) stop("not yet implemented")

        D1 <- L2Fam@distribution
        if(!is(D1, "DiscreteDistribution")) stop("not yet implemented")

        sign <- sign(biastype)
        w0 <- options("warn")
        on.exit(options(w0))
        options(warn = -1)
        L2deriv <- L2Fam@L2derivDistr[[1]]        
        
        l <- length(support(L2deriv))
        if (sign>0)
           {z0 <- support(L2deriv)[1]; deltahat <- support(L2deriv)[2]-z0
        }else{
            z0 <- support(L2deriv)[l]; deltahat <- z0-support(L2deriv)[l-1]}
        p0 <- d(L2deriv)(z0)   
        
        rad <- sqrt((abs(z0)/deltahat-(1-p0))/p0)
        names(rad) <- "lower case radius"

       options(w0)
       return(rad)
    })

# trick to make it callable from minmaxBias
setMethod("lowerCaseRadius", signature(L2Fam = "UnivariateDistribution",
                                       neighbor = "ContNeighborhood",
                                       risk = "asMSE",
                                       biastype = "onesidedBias"),
    function(L2Fam, neighbor, risk, biastype){

        L2deriv <- D1 <- L2Fam
        if(!is(D1, "DiscreteDistribution")) stop("not yet implemented")

        sign <- sign(biastype)
        w0 <- options("warn")
        on.exit(options(w0))
        options(warn = -1)
        
        l <- length(support(L2deriv))
        if (sign>0)
           {z0 <- support(L2deriv)[1]; deltahat <- support(L2deriv)[2]-z0
        }else{
            z0 <- support(L2deriv)[l]; deltahat <- z0-support(L2deriv)[l-1]}
        p0 <- d(L2deriv)(z0)   
        
        rad <- sqrt((abs(z0)/deltahat-(1-p0))/p0)
        names(rad) <- "lower case radius"

       options(w0)
       return(rad)
    })


 setMethod("lowerCaseRadius", signature(L2Fam = "L2ParamFamily",
                                       neighbor = "ContNeighborhood",
                                       risk = "asMSE",
                                       biastype = "asymmetricBias"),
    function(L2Fam, neighbor, risk, biastype){
        if(length(L2Fam@param) != 1) stop("not yet implemented")

        D1 <- L2Fam@distribution
        if(!is(D1, "DiscreteDistribution")) stop("not yet implemented")

        nu1 <- nu(biastype)[1]
        nu2 <- nu(biastype)[2]

        w0 <- options("warn")
        options(warn = -1)
        L2deriv <- L2Fam@L2derivDistr[[1]]        

        supp <- support(L2deriv)
        l <- length(supp)

        num <- nu2/(nu1+nu2)        
        
        zl <- q(L2deriv)(num)
        pl <- p(L2deriv)(zl)
        dl <- d(L2deriv)(zl)
        
        if (pl > num)
           { zm <- zu <- zl            
             wsm <- 0 
        
        } else {
            zu <- min(supp[p(L2deriv)[supp]>num])
             zm <- (zl*nu2+zu*nu1)/(nu1+nu2)
             wsm <- dl 
           }
        
        gg1 <- min(supp[supp > zm] - zm)
        gg2 <- max(supp[supp < zm] - zm)
        gg <- min(abs(supp[supp != zm] - zm))
    
        
        Int <- E(L2deriv, function(x, m){abs(x-m)}, m = zm)
        omega <- 2/(Int/nu1+Int/nu2)

        if(wsm > 0){
            if( (((zm == supp[1]) | (zm == supp[l])) & gg>0) |
                ((zm > supp[1]) & (zm < supp[l]) & (min(gg1,gg2)>0 ) ))
            {
            del1 <- pl-num
            del2 <- num-pl+dl
            M1 <- (del1*nu2*(nu1+1)+del2*nu1*(nu2-1))/
                      (del1+del2)/nu1^2/nu2/gg1
            M2 <- (del2*nu1*(nu2+1)+del2*nu2*(1-nu1))/
                      (del1+del2)/nu1/nu2^2/gg2
            M <- max(M1,M2)
            if (zm == supp[1]) M <- M1
            if (zm == supp[l]) M <- M2
            
            Int2 <- 1/nu1/nu2    

            }else{
                options(w0)
                rad <- Inf
                names(rad) <- "lower case radius"
                return(rad)
            }
        }else{
            M <- (nu1+nu2)/nu1/nu2/(zu-zl)
            ga <- ((pl-dl)/nu2-(1-pl)/nu1)/dl
            Int2 <- (1-pl)/nu1^2+(pl-dl)/nu2^2+dl*ga^2                                    
        }

       rad <- sqrt(M/omega- Int2)
       names(rad) <- "lower case radius"
       options(w0)
       return(rad)            
    })

