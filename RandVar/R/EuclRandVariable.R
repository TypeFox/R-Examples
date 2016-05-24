## generating function
EuclRandVariable <- function(Map = list(function(x){1}), 
                             Domain = NULL, dimension = 1, Range){
    nrvalues <- length(Map)
    for(i in 1:nrvalues){
        if(!is.function(Map[[i]])) 
            stop("element ", i, " of 'Map' contains no function")
        if(length(formals(Map[[i]])) != 1)
            stop("element ", i, " of 'Map' has to be a function of one argument")
        if(names(formals(Map[[i]])) != "x")
            stop("element ", i, " of 'Map' contains a function with argument name != 'x'")
    }

    R <- new("EuclRandVariable")
    R@Map <- Map
    R@Domain <- Domain

    if(missing(Range)) 
        R@Range <- EuclideanSpace(dimension = dimension)
    else{
        if(!is(Range, "EuclideanSpace"))
            stop("'Range' is no Euclidean space")
        R@Range <- Range
    }
    
    return(R)
}

## replace method
setReplaceMethod("Range", "EuclRandVariable", 
    function(object, value){ 
        object@Range <- value 
        if(!is(value, "EuclideanSpace"))
            stop("Range of 'value' is no Euclidean space")
        object
    })
    
## generating function
EuclRandMatrix <- function(Map = list(function(x){1}), nrow = 1, ncol = 1,
                              Domain = NULL, dimension = 1, Range) {
    nrvalues <- length(Map)
    for(i in 1:nrvalues){
        if(!is.function(Map[[i]])) 
            stop("element ", i, " of 'Map' contains no function")
        if(length(formals(Map[[i]])) != 1)
            stop("element ", i, " of 'Map' has to be a function of one argument")
        if(names(formals(Map[[i]])) != "x")
            stop("element ", i, " of 'Map' contains a function with argument name != 'x'")
    }
    if(missing(nrow)) 
        nrow <- ceiling(length(Map)/ncol)
    else if (missing(ncol)) 
        ncol <- ceiling(length(Map)/nrow)
    
    if(length(Map) != nrow*ncol)
        stop("'Map' has wrong dimension")

    R <- new("EuclRandMatrix")
    R@Map <- Map
    R@Domain <- Domain

    if(missing(Range)) 
        R@Range <- EuclideanSpace(dimension = dimension)
    else{
        if(!is(Range, "EuclideanSpace"))
            stop("'Range' is no Euclidean space")
        R@Range <- Range
    }
    
    R@Dim <- as.integer(c(nrow, ncol))
    
    return(R)
}

## access methods
setMethod("Dim", "EuclRandMatrix", function(object) object@Dim)
setMethod("nrow", "EuclRandMatrix", function(x) x@Dim[1])
setMethod("ncol", "EuclRandMatrix", function(x) x@Dim[2])

## setAs
setAs(from = "EuclRandVariable", to = "EuclRandMatrix", 
    def = function(from){
        R <- new("EuclRandMatrix") 
        R@Map <- from@Map
        R@Dim <- as.integer(c(length(from), 1))
        R@Domain <- from@Domain
        R@Range <- from@Range
        
        return(R)
    })

## replace methods
setReplaceMethod("Dim", "EuclRandMatrix", 
    function(object, value){ 
        d <- object@Dim
        val <- as.integer(value)
        if(d[1]*d[2] != val[1]*val[2])
            stop(paste("dims [ product", val[1]*val[2], "] do not match the length of object [", d[1]*d[2],"]"))
        object@Dim <- val
        object
    })

## dimension
setMethod("dimension", "EuclRandVariable", 
    function(object){ 
        length(object)*dimension(Range(object)) 
    })
setMethod("dimension", "EuclRandMatrix", 
    function(object){ 
        dimension(as(object, "EuclRandVariable"))
    })

## Extract via "["
setMethod("[", "EuclRandVariable",
    function(x, i, j, ..., drop = TRUE){
        if(!missing(j))
            stop("incorrect number of dimensions")
        x@Map <- x@Map[i]

        return(x)
    })
setMethod("[", "EuclRandMatrix",
    function(x, i, j, ..., drop = TRUE){
        map <- matrix(x@Map, nrow = x@Dim[1])
        if(missing(i))
            if(missing(j))
                map <- map[, ,..., drop = drop]
            else
                map <- map[, j, ..., drop = drop]
        else
            if(missing(j))
                map <- map[i, ,..., drop = drop]
            else
                map <- map[i, j, ..., drop = drop]
        if(is.matrix(map)){
            x@Map <- unlist(map)
            x@Dim <- as.integer(c(nrow(map), ncol(map)))
            return(x)
        }else{
            x@Map <- map            
            return(as(x, "EuclRandVariable"))
        }
    })


## evaluate Map of EuclRandVariable
setMethod("evalRandVar", signature(RandVar = "EuclRandVariable", 
                                   x = "numeric",
                                   distr = "missing"), 
    function(RandVar, x){
        if(!is.null(RandVar@Domain)){
            if(length(x) != RandVar@Domain@dimension)
                stop("x has wrong dimension")
        }else{
            warning("domain of 'RandVar' is 'NULL'")
        }

        nrvalues <- length(RandVar)
        res <- matrix(0, nrow = nrvalues, ncol = RandVar@Range@dimension)
        
        for(i in 1:nrvalues) res[i,] <- RandVar@Map[[i]](x)
        
        return(res)
    })
setMethod("evalRandVar", signature(RandVar = "EuclRandVariable", 
                                   x = "matrix",
                                   distr = "missing"), 
    function(RandVar, x){
        if(!is.null(RandVar@Domain)){
            if(ncol(x) != RandVar@Domain@dimension)
                stop("x has wrong dimension")
        }else{
            warning("domain of 'RandVar' is 'NULL'")
        }

        nrvalues <- length(RandVar)
        res <- array(0, c(nrvalues, nrow(x), RandVar@Range@dimension))
        
        for(i in 1:nrvalues){
            fun <- RandVar@Map[[i]]
            res[i,,] <- t(apply(x, 1, fun))
        }
        
        return(res)
    })
setMethod("evalRandVar", signature(RandVar = "EuclRandVariable", 
                                   x = "numeric",
                                   distr = "Distribution"), 
    function(RandVar, x, distr){
        if(is.null(RandVar@Domain))
            stop("domain of 'RandVar' is 'NULL'")
        if(length(x) != RandVar@Domain@dimension)
            stop("x has wrong dimension")
        if(RandVar@Domain@dimension != dimension(img(distr)))
            stop("'img' of 'distr' has wrong dimension")

        nrvalues <- length(RandVar)
        res <- matrix(NA, nrow = nrvalues, ncol = RandVar@Range@dimension)
        
        if(liesInSupport(distr, x))
            for(i in 1:nrvalues) res[i,] <- RandVar@Map[[i]](x)
        
        return(res)
    })
setMethod("evalRandVar", signature(RandVar = "EuclRandVariable", 
                                   x = "matrix",
                                   distr = "Distribution"), 
    function(RandVar, x, distr){
        if(is.null(RandVar@Domain))
            stop("domain of 'RandVar' is 'NULL'")
        if(ncol(x) != RandVar@Domain@dimension)
            stop("x has wrong dimension")
        if(RandVar@Domain@dimension != dimension(img(distr)))
            stop("'img' of 'distr' has wrong dimension")

        nrvalues <- length(RandVar)
        res <- array(NA, c(nrvalues, nrow(x), RandVar@Range@dimension))
        
        for(i in 1:nrvalues){
            fun <- RandVar@Map[[i]]
            for(j in 1:nrow(x))
                if(!liesInSupport(distr, x[j,]))
                    next
                else
                    res[i,j,] <- fun(x[j,])
        }
        
        return(res)
    })
setMethod("evalRandVar", signature(RandVar = "EuclRandMatrix", 
                                   x = "numeric",
                                   distr = "missing"), 
    function(RandVar, x){
        if(!is.null(RandVar@Domain)){
            if(length(x) != RandVar@Domain@dimension)
                stop("x has wrong dimension")
        }else{
            warning("domain of 'RandVar' is 'NULL'")
        }

        d <- RandVar@Dim
        res <- array(0, c(d[1], d[2], RandVar@Range@dimension))
        for(i in 1:d[1])
            for(j in 1:d[2])
                res[i,j,] <- RandVar@Map[[(i-1)*d[2] + j]](x)

        return(res)
    })
setMethod("evalRandVar", signature(RandVar = "EuclRandMatrix", 
                                   x = "matrix",
                                   distr = "missing"), 
    function(RandVar, x){
        if(!is.null(RandVar@Domain)){
            if(ncol(x) != RandVar@Domain@dimension)
                stop("x has wrong dimension")
        }else{
            warning("domain of 'RandVar' is 'NULL'")
        }

        d <- RandVar@Dim
        res <- array(0, c(d[1], d[2], nrow(x), RandVar@Range@dimension))
        
        for(i in 1:d[1]) 
            for(j in 1:d[2]){
                fun <- RandVar@Map[[(i-1)*d[2] + j]]
                res[i,j,,] <- t(apply(x, 1, fun))
            }
        
        return(res)
    })
setMethod("evalRandVar", signature(RandVar = "EuclRandMatrix", 
                                   x = "numeric",
                                   distr = "Distribution"), 
    function(RandVar, x, distr){
        if(is.null(RandVar@Domain))
            stop("domain of 'RandVar' is 'NULL'")
        if(ncol(x) != RandVar@Domain@dimension)
            stop("x has wrong dimension")
        if(RandVar@Domain@dimension != dimension(img(distr)))
            stop("'img' of 'distr' has wrong dimension")

        d <- RandVar@Dim
        res <- array(NA, c(d[1], d[2], RandVar@Range@dimension))

        if(liesInSupport(distr, x)){
            for(i in 1:d[1])
                for(j in 1:d[2])
                    res[i,j,] <- RandVar@Map[[(i-1)*d[2] + j]](x)
        }

        return(res)
    })
setMethod("evalRandVar", signature(RandVar = "EuclRandMatrix", 
                                   x = "matrix",
                                   distr = "Distribution"), 
    function(RandVar, x, distr){
        if(is.null(RandVar@Domain))
            stop("domain of 'RandVar' is 'NULL'")
        if(ncol(x) != RandVar@Domain@dimension)
            stop("x has wrong dimension")
        if(RandVar@Domain@dimension != dimension(img(distr)))
            stop("'img' of 'distr' has wrong dimension")

        d <- RandVar@Dim
        res <- array(NA, c(d[1], d[2], nrow(x), RandVar@Range@dimension))
        
        for(i in 1:d[1]) 
            for(j in 1:d[2]){
                fun <- RandVar@Map[[(i-1)*d[2] + j]]
                for(k in 1:nrow(x))
                    if(!liesInSupport(distr, x[k,]))
                        next
                    else
                        res[i,j,k,] <- fun(x[k,])
            }
        
        return(res)
    })

## computation of image distribution
setMethod("imageDistr", signature(RandVar = "EuclRandVariable", 
                                  distr = "Distribution"), 
    function(RandVar, distr){ 
        if(RandVar@Domain@dimension != dimension(img(distr)))
            stop("dimension of domain of 'RandVar' != dimension of img of 'distr'")

        nrvalues <- length(RandVar)
        res <- vector(mode = "list", length = nrvalues)
        for(i in 1:nrvalues)
            res[[i]] <- .getImageDistr(f = RandVar@Map[[i]], distr)
        return(new("DistrList", res))
    })

## transpose
setMethod("t", signature(x = "EuclRandVariable"), 
    function(x){ 
        nrvalues <- length(x)
        map <- vector(mode = "list", length = nrvalues)
        fct <- NULL
        for(i in 1:nrvalues){
            map[[i]] <- function(x){ f <- fct; t(f(x)) }
            body(map[[i]]) <- substitute({ f <- fct; t(f(x)) },
                                    list(fct = x@Map[[i]]))
        }

        R <- new("EuclRandMatrix") 
        R@Map <- map
        R@Dim <- as.integer(c(1, length(x)))
        R@Domain <- x@Domain
        R@Range <- x@Range
        
        return(R)
    })
setMethod("t", signature(x = "EuclRandMatrix"), 
    function(x){ 
        map <- matrix(x@Map, nrow = x@Dim[1])

        fkt <- NULL
        d <- x@Dim
        map <- vector(mode = "list", length = d[1]*d[2])
        for(i in 1:d[1])
            for(j in 1:d[2]){
                map[[(i-1)*d[2] + j]] <- function(x){ f <- fkt; t(f(x)) }
                body(map[[(i-1)*d[2] + j]]) <- substitute({ f <- fkt; t(f(x)) },
                                                    list(fkt = x@Map[[(j-1)*d[1] + i]]))
        }
        x@Map <- map
        x@Dim <- as.integer(c(d[2], d[1]))

        return(x)
    })
