# generating function for class 'EuclRandVarList'
EuclRandVarList <- function(...){ 
    new("EuclRandVarList", list(...)) 
}

setAs(from = "EuclRandVariable", to = "EuclRandVarList", 
    def = function(from){ new("EuclRandVarList", list(from)) })

setMethod("dimension", "EuclRandVarList", 
    function(object){ 
        dimn <- 0
        for(i in 1:length(object))
            dimn <- dimn + dimension(object[[i]])
        return(dimn)
    })
setMethod("numberOfMaps", "EuclRandVarList", 
    function(object){ 
        nr <- 0
        for(i in 1:length(object))
            nr <- nr + length(object[[i]])
        return(nr)
    })

setMethod("t", "EuclRandVarList",
    function(x){ for(i in 1:length(x)) x[[i]] <- t(x[[i]]) })

setMethod("evalRandVar", signature(RandVar = "EuclRandVarList", 
                                   x = "numeric",
                                   distr = "missing"), 
    function(RandVar, x){
        if(!is.null(RandVar[[1]]@Domain)){
            if(length(x) != RandVar[[1]]@Domain@dimension)
                stop("x has wrong dimension")
        }else{
            warning("domain of 'RandVar' is 'NULL'")
        }

        nrvalues <- length(RandVar)
        res <- vector("list", nrvalues)
        
        for(i in 1:nrvalues) res[[i]] <- evalRandVar(RandVar[[i]], x)
        
        return(res)
    })
setMethod("evalRandVar", signature(RandVar = "EuclRandVarList", 
                                   x = "matrix",
                                   distr = "missing"), 
    function(RandVar, x){
        if(!is.null(RandVar[[1]]@Domain)){
            if(ncol(x) != RandVar[[1]]@Domain@dimension)
                stop("x has wrong dimension")
        }else{
            warning("domain of 'RandVar' is 'NULL'")
        }

        nrvalues <- length(RandVar)
        res <- vector("list", nrvalues)
        
        for(i in 1:nrvalues) res[[i]] <- evalRandVar(RandVar[[i]], x)
        
        return(res)
    })
setMethod("evalRandVar", signature(RandVar = "EuclRandVarList", 
                                   x = "numeric",
                                   distr = "Distribution"), 
    function(RandVar, x, distr){
        if(is.null(RandVar[[1]]@Domain))
            stop("domain of 'RandVar' is 'NULL'")
        if(ncol(x) != RandVar[[1]]@Domain@dimension)
            stop("x has wrong dimension")
        if(RandVar[[1]]@Domain@dimension != dimension(img(distr)))
            stop("x has wrong dimension != dimension of 'img' of 'distr'")

        nrvalues <- length(RandVar)
        res <- vector("list", nrvalues)
        
        for(i in 1:nrvalues) res[[i]] <- evalRandVar(RandVar[[i]], x, distr)
        
        return(res)
    })
setMethod("evalRandVar", signature(RandVar = "EuclRandVarList", 
                                   x = "matrix",
                                   distr = "Distribution"), 
    function(RandVar, x, distr){
        if(is.null(RandVar[[1]]@Domain))
            stop("domain of 'RandVar' is 'NULL'")
        if(ncol(x) != RandVar[[1]]@Domain@dimension)
            stop("x has wrong dimension")
        if(RandVar[[1]]@Domain@dimension != dimension(img(distr)))
            stop("x has wrong dimension != dimension of 'img' of 'distr'")

        nrvalues <- length(RandVar)
        res <- vector("list", nrvalues)
        
        for(i in 1:nrvalues) res[[i]] <- evalRandVar(RandVar[[i]], x, distr)
        
        return(res)
    })

setMethod("imageDistr", signature(RandVar = "EuclRandVarList", 
                                  distr = "Distribution"), 
    function(RandVar, distr){ 
        if(RandVar[[1]]@Domain@dimension != dimension(img(distr)))
            stop("dimension of domain of 'RandVar' != dimension of 'img' of 'distr'")

        nrvalues1 <- length(RandVar)
        nrvalues2 <- 0
        for(i in 1:length(nrvalues1))
            nrvalues2 <- nrvalues2 + length(RandVar[[i]])
        
        res <- vector("list", nrvalues2)
        comp <- 0
        for(i in 1:nrvalues1){
            for(j in 1:length(RandVar[[i]])){
                comp <- comp + 1
                res[[comp]] <-  .getImageDistr(f = RandVar[[i]]@Map[[j]], distr)
            }
        }

        return(new("DistrList", res))
    })

## matrix multiplication
setMethod("%m%", signature(x = "EuclRandVarList", y = "EuclRandVarList"),
    function(x, y){
        nrvalues <- length(x)
        if(nrvalues != length(y))
            stop("non-conformable arguments")
        
        res <- vector("list", nrvalues^2)
        comp <- 0
        for(i in 1:nrvalues)
            for(j in 1:nrvalues){
                comp <- comp + 1
                res[[comp]] <- x[[j]] %*% t(y[[i]])
        }
        x@.Data <- res

        return(x)
    })
