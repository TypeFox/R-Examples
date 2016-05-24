## generate IC
## for internal use only!
setMethod("generateIC", signature(neighbor = "TotalVarNeighborhood", 
                                  L2Fam = "L2RegTypeFamily"),
    function(neighbor, L2Fam, res){
        A <- res$A
        a <- sign(as.vector(A))*res$a
        b <- res$b
        ICfct <- vector(mode = "list", length = 1)
        Y <- as(A %*% L2Fam@L2deriv, "EuclRandVariable")
        if(!is.null(res$d)){
            a <- as.vector(A)*a
            ICfct[[1]] <- function(x){ ind1 <- (Y(x) > 0); ind2 <- (Y(x) < 0)
                                       (a+b)*ind1 + a*ind2 }
            body(ICfct[[1]]) <- substitute({ ind1 <- (Y(x) > 0); ind2 <- (Y(x) < 0)
                                             (a+b)*ind1 + a*ind2 },
                                             list(Y = Y@Map[[1]], a = a, b = b))
        }else{
            if((a == -Inf) & (b == Inf)){
                ICfct[[1]]<- function(x){ Y(x) }
                body(ICfct[[1]]) <- substitute({ Y(x) }, list(Y = Y@Map[[1]]))
            }else{
                ICfct[[1]] <- function(x){ pmin(pmax(a, Y(x)), a+b) }
                body(ICfct[[1]]) <- substitute({ pmin(pmax(a, Y(x)), a+b) },
                                                 list(Y = Y@Map[[1]], a = a, b = b))
            }
        }
        if((a == -Inf) & (b == Inf))
            clipUp <- Inf
        else
            clipUp <- a + b
        return(TotalVarIC(
                name = "IC of total variation type", 
                CallL2Fam = call("L2RegTypeFamily", 
                                name = L2Fam@name,
                                distribution = L2Fam@distribution,  
                                param = L2Fam@param,
                                props = L2Fam@props,
                                ErrorDistr = L2Fam@ErrorDistr,
                                ErrorSymm = L2Fam@ErrorSymm,
                                RegDistr = L2Fam@RegDistr,
                                RegSymm = L2Fam@RegSymm,
                                Regressor = L2Fam@Regressor,
                                L2deriv = L2Fam@L2deriv,
                                ErrorL2deriv = L2Fam@ErrorL2deriv,
                                ErrorL2derivDistr = L2Fam@ErrorL2derivDistr,
                                ErrorL2derivSymm = L2Fam@ErrorL2derivSymm,
                                FisherInfo = L2Fam@FisherInfo),
                Curve = EuclRandVarList(EuclRandVariable(Map = ICfct, Domain = Y@Domain, 
                                         Range = Y@Range)),
                clipUp = clipUp,
                clipLo = a,
                stand = A,
                neighborRadius = neighbor@radius,
                Risks = res$risk,
                Infos = matrix(res$info, ncol = 2, 
                            dimnames = list(character(0), c("method", "message")))))
    })
