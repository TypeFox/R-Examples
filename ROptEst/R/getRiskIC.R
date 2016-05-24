###############################################################################
## asymptotic covariance
###############################################################################
setMethod("getRiskIC", signature(IC = "HampIC", 
                                 risk = "asCov",
                                 neighbor = "missing",
                                 L2Fam = "missing"),
    function(IC, risk){
        L2Fam <- force(eval(IC@CallL2Fam))
        getRiskIC(IC = IC, risk = risk, L2Fam = L2Fam)
    })

setMethod("getRiskIC", signature(IC = "HampIC", 
                                 risk = "asCov",
                                 neighbor = "missing",
                                 L2Fam = "L2ParamFamily"),
    function(IC, risk, L2Fam){
        Cov <- eval(IC@Risks[["asCov"]])
        if(is.null(Cov)){
           if(numberOfMaps(L2Fam@L2deriv)==1){ ## L2derivDim <- L2Fam@L2deriv
              L2deriv <- L2Fam@L2derivDistr[[1]]
              A <- as.vector(IC@stand)
              c0 <- IC@clip/abs(A)
              z <- IC@cent/A
              neighbor <- ContNeighborhood(1)
              Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor,
                         biastype = biastype(IC), clip = c0, cent = z, stand = A)
              return(list(asCov = list(distribution = .getDistr(L2Fam),
                          value = Cov)))
            }
            return(getRiskIC(as(IC, "IC"), risk = risk, L2Fam = L2Fam))
        }else
            return(list(asCov = list(distribution = .getDistr(L2Fam), value = Cov)))
    })

setMethod("getRiskIC", signature(IC = "TotalVarIC", 
                                 risk = "asCov",
                                 neighbor = "missing",
                                 L2Fam = "L2ParamFamily"),
    function(IC, risk, L2Fam){
        Cov <- eval(IC@Risks[["asCov"]])
        if (is.null(Cov)){
            L2deriv <- L2Fam@L2derivDistr[[1]]
            A <- as.vector(IC@stand)
            c0 <- (IC@clipUp-IC@clipLo)/abs(A)
            z <- IC@clipLo/abs(A)
            neighbor <- TotalVarNeighborhood(1)
            Cov <- getInfV(L2deriv = L2deriv, neighbor = neighbor, 
                       biastype = biastype(IC), clip = c0, cent = z, stand = A)
            }
        return(list(asCov = list(distribution = .getDistr(L2Fam), value = Cov)))
    })

###############################################################################
## asymptotic Bias for various types
###############################################################################
setMethod("getBiasIC", signature(IC = "HampIC",
                                 neighbor = "UncondNeighborhood"),
    function(IC, neighbor, L2Fam,...){
        if(missing(L2Fam))
            L2Fam <- force(eval(IC@CallL2Fam))

            Bias <- IC@Risks$asBias$value

        return(list(asBias = list(distribution = .getDistr(L2Fam), 
                    neighborhood = neighbor@type, value = Bias)))
    })

setMethod("getBiasIC", signature(IC = "TotalVarIC",
                                 neighbor = "UncondNeighborhood"),
    function(IC, neighbor, L2Fam,...){
        if(missing(L2Fam))
            L2Fam <- force(eval(IC@CallL2Fam))

        Bias <- IC@Risks$asBias$value
        if (is.null(Bias)){
            Bias <- if(is(neighbor,"ContNeighborhood"))
                     max(IC@clipUp,-IC@clipLo) else IC@clipUp-IC@clipLo
        }

        return(list(asBias = list(distribution = .getDistr(L2Fam), 
                    neighborhood = neighbor@type, value = Bias)))
    })

