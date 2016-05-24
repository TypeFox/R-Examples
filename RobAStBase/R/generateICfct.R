## generate IC
## for internal use only!
setMethod("generateIC.fct", signature(neighbor = "UncondNeighborhood", L2Fam = "L2ParamFamily"),
    function(neighbor, L2Fam, res){
        A <- as.matrix(res$A)
        a <- if(is(neighbor,"TotalVarNeighborhood")) 0 else res$a 
        b <- res$b
        d <- res$d
        w <- weight(res$w)
        nrvalues <- nrow(A)
        dim <- ncol(A)
        ICfct <- vector(mode = "list", length = nrvalues)
        Y <- as(A %*% L2Fam@L2deriv - a, "EuclRandVariable")
        L <- as(diag(dim)%*%L2Fam@L2deriv, "EuclRandVariable")
        L.fct <- function(x) evalRandVar(L,x)
        if(nrvalues == 1){
            if(!is.null(d)){
                ICfct[[1]] <- function(x){}
                if(all(dim(trafo(L2Fam@param)) == c(1, 1))){
                    body(ICfct[[1]]) <- substitute(
                                            { ind <- 1-.eq(Y(x))
                                              Y(x)*w(L(x)) + zi*(1-ind)*d*b },
                                            list(Y = Y@Map[[1]], L = L.fct, w = w, b = b, d = d,
                                                zi = sign(trafo(L2Fam@param)), .eq = .eq))
                }else{
                    body(ICfct[[1]]) <- substitute(
                                            { ind <- 1-.eq(Y(x))
                                              ifelse(ind, Y(x)*w(L(x)), NA) },
                                            list(Y = Y@Map[[1]], L = L.fct, w = w, b = b, d = d, 
                                                 .eq = .eq))
                }
            }else{
                ICfct[[1]] <- function(x){}
                body(ICfct[[1]]) <- substitute({ Y(x)*w(L(x)) },
                                                 list(Y = Y@Map[[1]], L = L.fct, w = w))
            }
        }else{
            if(!is.null(d))
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){}
                    body(ICfct[[i]]) <- substitute({ind <- 1-.eq(Yi(x))
                                                    ind*Yi(x)*w(L(x)) + (1-ind)*d
                                                    },
                                                 list(Yi = Y@Map[[i]], L = L.fct, w = w,
                                                      b = b, d = d[i]))#,  .eq = .eq))
                }
            else
                for(i in 1:nrvalues){
                    ICfct[[i]] <- function(x){}
                    body(ICfct[[i]]) <- substitute({  Yi(x)*w(L(x))  },
                                                 list(Yi = Y@Map[[i]], L = L.fct, w = w))
                }
        }
        return(EuclRandVarList(EuclRandVariable(Map = ICfct, Domain = Y@Domain,
                                         Range = Y@Range)))
    })

