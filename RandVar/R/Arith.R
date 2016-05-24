## "Arith" group
setMethod("Arith", signature(e1 = "numeric", e2 = "EuclRandVariable"),
    function(e1, e2){
        nrvalues1 <- length(e1)
        dimn <- e2@Range@dimension
        nrvalues2 <- length(e2)*dimn
        if(nrvalues1 != nrvalues2){
            if(nrvalues2 > nrvalues1){
                if(nrvalues2 %% nrvalues1 == 0){
                    e1 <- rep(e1, nrvalues2/nrvalues1)
                }else{
                    e1 <- rep(e1, nrvalues2/nrvalues1 + 1)
                    warning("longer object length is not a multiple of shorter object length")
                }
            }else{
                stop("length of 'numeric' has to be less or equal dimension of 'EuclRandVariable'")
            }
        }

        fct <- NULL;  f <- function(x,y){}
        map <- vector("list", length(e2))
        for(i in 1:length(e2)){
            map[[i]] <- function(x){ f2 <- fct; f(e1, f2(x)) }
            body(map[[i]]) <- substitute({ f2 <- fct; f(e1, f2(x)) },
                                    list(f = as.name(.Generic), 
                                         fct = e2@Map[[i]], e1 = e1[((i-1)*dimn + 1):(i*dimn)]))
        }
        e2@Map <- map

        return(e2)
    })
setMethod("Arith", signature(e1 = "numeric", e2 = "EuclRandMatrix"),
    function(e1, e2){
        e2@Map <- Map(callGeneric(e1, as(e2, "EuclRandVariable")))

        return(e2)
    })
setMethod("Arith", signature(e1 = "numeric", e2 = "EuclRandVarList"),
    function(e1, e2){
        nrvalues1 <- length(e1)
        nrvalues2 <- numeric(length(e2))
        for(i in 1:length(e2))
            nrvalues2[i] <- e2[[i]]@Range@dimension*length(e2[[i]])
        
        nrvalues3 <- sum(nrvalues2)
        
        if(nrvalues1 != nrvalues3){
            if(nrvalues3 > nrvalues1){
                if(nrvalues3 %% nrvalues1 == 0){
                    e1 <- rep(e1, nrvalues3/nrvalues1)
                }else{
                    e1 <- rep(e1, nrvalues3/nrvalues1 + 1)
                    warning("longer object length is not a multiple of shorter object length")
                }
            }else{
                stop("length of 'numeric' has to be less or equal dimension of 'EuclRandVarList'")
            }
        }

        for(i in 1:length(e2))
            e2[[i]] <- callGeneric(e1[((i-1)*nrvalues2[i] + 1):(i*nrvalues2[i])], e2[[i]])
        
        return(e2)
    })
setMethod("Arith", signature(e1 = "EuclRandVariable", e2 = "numeric"),
    function(e1, e2){
        dimn <- e1@Range@dimension
        nrvalues1 <- length(e1)*dimn
        nrvalues2 <- length(e2)
        if(nrvalues1 != nrvalues2){
            if(nrvalues1 > nrvalues2){
                if(nrvalues1 %% nrvalues2 == 0){
                    e2 <- rep(e2, nrvalues1/nrvalues2)
                }else{
                    e2 <- rep(e2, nrvalues1/nrvalues2 + 1)
                    warning("longer object length is not a multiple of shorter object length")
                }
            }else{
                stop("length of 'numeric' has to be less or equal dimension of 'EuclRandVariable'")
            }
        }

        fct <- NULL;  f <- function(x,y){}
        map <- vector("list", length(e1))
        for(i in 1:length(e1)){
            map[[i]] <- function(x){ f1 <- fct; f(f1(x), e2) }
            body(map[[i]]) <- substitute({ f1 <- fct; f(f1(x), e2) },
                                    list(f = as.name(.Generic), 
                                         fct = e1@Map[[i]], e2 = e2[((i-1)*dimn + 1):(i*dimn)]))
        }
        e1@Map <- map

        return(e1)
    })
setMethod("Arith", signature(e1 = "EuclRandMatrix", e2 = "numeric"),
    function(e1, e2){
        e1@Map <- Map(callGeneric(as(e1, "EuclRandVariable"), e2))

        return(e1)
    })
setMethod("Arith", signature(e1 = "EuclRandVarList", e2 = "numeric"),
    function(e1, e2){
        nrvalues1 <- numeric(length(e1))
        for(i in 1:length(e1))
            nrvalues1[i] <- e1[[i]]@Range@dimension*length(e1[[i]])

        nrvalues2 <- length(e2)        
        nrvalues3 <- sum(nrvalues1)
        
        if(nrvalues2 != nrvalues3){
            if(nrvalues3 > nrvalues2){
                if(nrvalues3 %% nrvalues2 == 0){
                    e2 <- rep(e2, nrvalues3/nrvalues2)
                }else{
                    e2 <- rep(e2, nrvalues3/nrvalues2 + 1)
                    warning("longer object length is not a multiple of shorter object length")
                }
            }else{
                stop("length of 'numeric' has to be less or equal dimension of 'EuclRandVarList'")
            }
        }

        for(i in 1:length(e1))
            e1[[i]] <- callGeneric(e1[[i]], e2[((i-1)*nrvalues1[i] + 1):(i*nrvalues1[i])])
        
        return(e1)
    })
setMethod("Arith", signature(e1 = "EuclRandVariable", e2 = "EuclRandVariable"),
    function(e1, e2){
        if(!compatibleDomains(e1, e2))
            stop("the domains are not compatible")

        nrvalues1 <- length(e1)
        nrvalues2 <- length(e2)
        if(nrvalues1 != nrvalues2){
            if(nrvalues1 > nrvalues2){
                if(nrvalues1 %% nrvalues2 == 0){
                    e2@Map <- rep(e2@Map, nrvalues1/nrvalues2)
                }else{
                    e2@Map <- rep(e2@Map, nrvalues1/nrvalues2 + 1)
                    warning("longer object length is not a multiple of shorter object length")
                }
            }else{
                if(nrvalues2 %% nrvalues1 == 0){
                    e1@Map <- rep(e1@Map, nrvalues2/nrvalues1)
                }else{
                    e1@Map <- rep(e1@Map, nrvalues2/nrvalues1 + 1)
                    warning("longer object length is not a multiple of shorter object length")
                }
            }
        }

        fct1 <- NULL;  fct2 <- NULL; f <- function(x,y){}
        nrvalues <- max(nrvalues1, nrvalues2)
        map <- vector("list", nrvalues)
        for(i in 1:nrvalues){
            map[[i]] <- function(x){ f1 <- fct1; f2 <- fct2; f(f1(x), f2(x)) }
            body(map[[i]]) <- substitute({ f1 <- fct1; f2 <- fct2; f(f1(x), f2(x)) },
                                    list(f = as.name(.Generic), 
                                         fct1 = e1@Map[[i]], fct2 = e2@Map[[i]]))
        }

        if(identical(e1@Domain, e2@Domain))
            Dom <- e1@Domain
        else
            if(is(e1@Domain, class(e2@Domain)[1]))
                Dom <- e1@Domain
            else
                Dom <- e2@Domain

        if(identical(e1@Range, e2@Range))
            Ran <- e1@Range
        else
            Ran <- EuclideanSpace(dimension = max(e1@Range@dimension, e2@Range@dimension))
            
        e1@Map <- map
        e1@Domain <- Dom
        e1@Range <- Ran

        return(e1)
    })
setMethod("Arith", signature(e1 = "EuclRandMatrix", e2 = "EuclRandMatrix"),
    function(e1, e2){
        if(!identical(e1@Dim, e2@Dim))
            stop("non-comformable random matrices")
        res <- callGeneric(as(e1, "EuclRandVariable"), as(e2, "EuclRandVariable"))
        e1@Map <- res@Map
        e1@Domain <- res@Domain
        e1@Range <- res@Range

        return(e1)
    })
setMethod("Arith", signature(e1 = "EuclRandVarList", e2 = "EuclRandVarList"),
    function(e1, e2){
        if(length(e1) != length(e2))
            stop("non-comformable lists of random variables")
        for(i in 1:length(e1))
            e1[[i]] <- callGeneric(e1[[i]], e2[[i]])
        
        return(e1)
    })
