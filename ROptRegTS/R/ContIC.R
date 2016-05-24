## generate IC
## for internal use only!
setMethod("generateIC", signature(neighbor = "ContNeighborhood", 
                                  L2Fam = "L2RegTypeFamily"),
    function(neighbor, L2Fam, res){
        A <- res$A
        a <- res$a
        b <- res$b
        d <- res$d
        nrvalues <- nrow(A)
        ICfct <- vector(mode = "list", length = nrvalues)
        Y <- as(A %*% L2Fam@L2deriv - a, "EuclRandVariable")
        if(nrvalues == 1){
            if(!is.null(d)){
                ICfct[[1]] <- function(x){ 
                                    ind <- (Y(x) != 0) 
                                    b*(ind*Y(x)/(ind*absY(x) + (1-ind)*1) + zi*(1-ind)*d)
                              }
                body(ICfct[[1]]) <- substitute(
                                        { ind <- (Y(x) != 0) 
                                          b*(ind*Y(x)/(ind*absY(x) + (1-ind)*1) + zi*(1-ind)*d) },
                                        list(Y = Y@Map[[1]], absY = abs(Y)@Map[[1]], b = b, d = d, 
                                             zi = sign(L2Fam@param@trafo)))
            }else{
                ICfct[[1]] <- function(x){ Y(x)*pmin(1, b/absY(x)) }
                body(ICfct[[1]]) <- substitute({ Y(x)*pmin(1, b/absY(x)) },
                                                 list(Y = Y@Map[[1]], absY = abs(Y)@Map[[1]], b = b))
            }
        }
        else{
            absY <- sqrt(Y %*% Y)
            if(!is.null(d))
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){ ind <- (Yi(x) != 0) ; ind*b*Yi(x)/absY(x) + (1-ind)*d }
                    body(ICfct[[i]]) <- substitute({ ind <- (Yi(x) != 0) ; ind*b*Yi(x)/absY(x) + (1-ind)*d },
                                                 list(Yi = Y@Map[[i]], absY = absY@Map[[1]], b = b, d = d[i]))
                }
            else
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){ Yi(x)*pmin(1, b/absY(x)) }
                    body(ICfct[[i]]) <- substitute({ Yi(x)*pmin(1, b/absY(x)) },
                                                 list(Yi = Y@Map[[i]], absY = absY@Map[[1]], b = b))
                }
        }
        return(ContIC(
               name = "IC of contamination type", 
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
                clip = b,
                cent = a,
                stand = A,
                lowerCase = d,
                neighborRadius = neighbor@radius,
                Risks = res$risk,
                Infos = matrix(res$info, ncol = 2, 
                            dimnames = list(character(0), c("method", "message")))))
    })
