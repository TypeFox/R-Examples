setMethod("E", signature(object = "UnivariateDistribution", 
                         fun = "EuclRandVariable", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        if(!is(fun@Domain, "EuclideanSpace"))
            stop("'Domain' of the random variable is no Euclidean space")
        if(dimension(fun@Domain) != 1)
            stop("dimension of 'Domain' of the random variable has to be 1")
        if(dimension(fun@Range) != 1)
            stop("dimension of 'Range' of the random variable has to be 1")

        dimn <- length(fun)
        res <- numeric(dimn)
        for(i in 1:dimn) res[i] <- E(object, fun = Map(fun)[[i]], useApply = useApply, ...)

        return(res)
    })
setMethod("E", signature(object = "AbscontDistribution", 
                         fun = "EuclRandVariable", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        if(!is(fun@Domain, "EuclideanSpace"))
            stop("'Domain' of the random variable is no Euclidean space")
        if(dimension(fun@Domain) != 1)
            stop("dimension of 'Domain' of the random variable has to be 1")
        if(dimension(fun@Range) != 1)
            stop("dimension of 'Range' of the random variable has to be 1")

        dimn <- length(fun)
        res <- numeric(dimn)
        for(i in 1:dimn) res[i] <- E(object, fun = Map(fun)[[i]], useApply = useApply, ...)

        return(res)
    })
setMethod("E", signature(object = "DiscreteDistribution", 
                         fun = "EuclRandVariable", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        if(!is(fun@Domain, "EuclideanSpace"))
            stop("'Domain' of the random variable is no Euclidean space")
        if(dimension(fun@Domain) != 1)
            stop("dimension of 'Domain' of the random variable has to be 1")
        if(dimension(fun@Range) != 1)
            stop("dimension of 'Range' of the random variable has to be 1")

        dimn <- length(fun)
        res <- numeric(dimn)
        for(i in 1:dimn) res[i] <- E(object, fun = Map(fun)[[i]], useApply = useApply, ...)

        return(res)
    })
setMethod("E", signature(object = "UnivariateDistribution", 
                         fun = "EuclRandMatrix", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        matrix(E(object, as(fun, "EuclRandVariable"), useApply = useApply, ...), nrow = nrow(fun))
    })
setMethod("E", signature(object = "AbscontDistribution", 
                         fun = "EuclRandMatrix", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        matrix(E(object, as(fun, "EuclRandVariable"), useApply = useApply, ...), nrow = nrow(fun))
    })
setMethod("E", signature(object = "DiscreteDistribution", 
                         fun = "EuclRandMatrix", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        matrix(E(object, as(fun, "EuclRandVariable"), useApply = useApply, ...), nrow = nrow(fun))
    })
setMethod("E", signature(object = "UnivariateDistribution", 
                         fun = "EuclRandVarList", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        nrvalues <- length(fun)
        res <- vector("list", nrvalues)
        for(i in 1:nrvalues) res[[i]] <- E(object, fun[[i]], useApply = useApply, ...)

        return(res)
    })
setMethod("E", signature(object = "AbscontDistribution", 
                         fun = "EuclRandVarList", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        nrvalues <- length(fun)
        res <- vector("list", nrvalues)
        for(i in 1:nrvalues) res[[i]] <- E(object, fun[[i]], useApply = useApply, ...)

        return(res)
    })
setMethod("E", signature(object = "DiscreteDistribution", 
                         fun = "EuclRandVarList", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        nrvalues <- length(fun)
        res <- vector("list", nrvalues)
        for(i in 1:nrvalues) res[[i]] <- E(object, fun[[i]], useApply = useApply, ...)

        return(res)
    })
setMethod("E", signature(object = "MultivariateDistribution", 
                         fun = "EuclRandVariable", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        if(!is(fun@Domain, "EuclideanSpace"))
            stop("'Domain' of the random variable is no Euclidean space")
        if(fun@Domain@dimension != object@img@dimension)
            stop("dimension of 'Domain' of the random variable is not equal\n",
                 "to dimension of 'img' of the distribution")
        dimn <- length(fun)
        res <- matrix(0, nrow = dimn, ncol = fun@Range@dimension)
        for(i in 1:dimn) res[i,] <- E(object, fun@Map[[i]], useApply = useApply, ...)

        return(res)
    })
setMethod("E", signature(object = "DiscreteMVDistribution", 
                         fun = "EuclRandVariable", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        if(!is(fun@Domain, "EuclideanSpace"))
            stop("'Domain' of the random variable is no Euclidean space")
        if(fun@Domain@dimension != object@img@dimension)
            stop("dimension of 'Domain' of the random variable is not equal\n",
                 "to dimension of 'img' of the distribution")
        dimn <- length(fun)
        res <- matrix(0, nrow = dimn, ncol = fun@Range@dimension)
        for(i in 1:dimn) res[i,] <- E(object, fun@Map[[i]], useApply = useApply, ...)

        return(res)
    })
setMethod("E", signature(object = "MultivariateDistribution", 
                         fun = "EuclRandMatrix", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        array(E(object, as(fun, "EuclRandVariable"), useApply = useApply, ...),
              c(nrow(fun), ncol(fun), fun@Range@dimension))
    })
setMethod("E", signature(object = "DiscreteMVDistribution", 
                         fun = "EuclRandMatrix", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        array(E(object, as(fun, "EuclRandVariable"), useApply = useApply, ...),
              c(nrow(fun), ncol(fun), fun@Range@dimension))
    })
setMethod("E", signature(object = "MultivariateDistribution", 
                         fun = "EuclRandVarList", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        nrvalues <- length(fun)
        res <- vector("list", nrvalues)
        for(i in 1:nrvalues)
            res[[i]] <- E(object, fun[[i]], useApply = useApply, ...)

        return(res)
    })
setMethod("E", signature(object = "DiscreteMVDistribution", 
                         fun = "EuclRandVarList", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE, ...){
        nrvalues <- length(fun)
        res <- vector("list", nrvalues)
        for(i in 1:nrvalues)
            res[[i]] <- E(object, fun[[i]], useApply = useApply, ...)

        return(res)
    })
setMethod("E", signature(object = "UnivariateCondDistribution", 
                         fun = "EuclRandVariable", 
                         cond = "numeric"),
    function(object, fun, cond, withCond = FALSE, useApply = TRUE, ...){
        if(!is(fun@Domain, "EuclideanSpace"))
            stop("'Domain' of the random variable has to be a Euclidean Space")
        if(withCond){
            if(fun@Domain@dimension != (1+length(cond)))
                stop("wrong dimension of 'Domain' of 'fun'")
        }else{
            if(fun@Domain@dimension != 1)
                stop("dimension of 'Domain' of 'fun' has to be 1")
        }
        if(dimension(fun@Range) != 1)
            stop("dimension of 'Range' of the random variable has to be 1")

        dimn <- length(fun)
        res <- numeric(dimn)
        if(withCond){
            for(i in 1:dimn){ 
                fun1 <- function(x, cond, fct){ fct(c(x, cond)) }
                res[i] <- E(object, fun1, cond, fct = fun@Map[[i]], 
                            withCond, useApply = useApply, ...)
            }
        }else{
            for(i in 1:dimn) res[i] <- E(object, fun@Map[[i]], cond, useApply = useApply, ...)
        }

        return(res)
    })
setMethod("E", signature(object = "AbscontCondDistribution", 
                         fun = "EuclRandVariable", 
                         cond = "numeric"),
    function(object, fun, cond, withCond = FALSE, useApply = TRUE, ...){
        if(!is(fun@Domain, "EuclideanSpace"))
            stop("'Domain' of the random variable has to be a Euclidean Space")
        if(withCond){
            if(fun@Domain@dimension != (1+length(cond)))
                stop("wrong dimension of 'Domain' of 'fun'")
        }else{
            if(fun@Domain@dimension != 1)
                stop("dimension of 'Domain' of 'fun' has to be 1")
        }
        if(dimension(fun@Range) != 1)
            stop("dimension of 'Range' of the random variable has to be 1")

        dimn <- length(fun)
        res <- numeric(dimn)
        if(withCond){
            for(i in 1:dimn){ 
                fun1 <- function(x, cond, fct){ fct(c(x, cond)) }
                res[i] <- E(object, fun1, cond, fct = fun@Map[[i]], 
                            withCond, useApply = useApply, ...)
            }
        }else{
            for(i in 1:dimn) res[i] <- E(object, fun@Map[[i]], cond, useApply = useApply, ...)
        }

        return(res)
    })
setMethod("E", signature(object = "DiscreteCondDistribution", 
                         fun = "EuclRandVariable", 
                         cond = "numeric"),
    function(object, fun, cond, withCond = FALSE, useApply = TRUE, ...){
        if(!is(fun@Domain, "EuclideanSpace"))
            stop("'Domain' of the random variable has to be a Euclidean Space")
        if(withCond){
            if(fun@Domain@dimension != (1+length(cond)))
                stop("wrong dimension of 'Domain' of 'fun'")
        }else{
            if(fun@Domain@dimension != 1)
                stop("dimension of 'Domain' of 'fun' has to be 1")
        }
        if(dimension(fun@Range) != 1)
            stop("dimension of 'Range' of the random variable has to be 1")

        dimn <- length(fun)
        res <- numeric(dimn)
        if(withCond){
            for(i in 1:dimn){ 
                fun1 <- function(x, cond, fct){ fct(c(x, cond)) }
                res[i] <- E(object, fun1, cond, fct = fun@Map[[i]], 
                            withCond, useApply = useApply, ...)
            }
        }else{                                                              
            for(i in 1:dimn) res[i] <- E(object, fun@Map[[i]], cond, useApply = useApply, ...)
        }

        return(res)
    })
setMethod("E", signature(object = "UnivariateCondDistribution",
                         fun = "EuclRandMatrix", 
                         cond = "numeric"),
    function(object, fun, cond, withCond = FALSE, useApply = TRUE, ...){
        matrix(E(object, as(fun, "EuclRandVariable"), cond, withCond, 
                 useApply = useApply, ...), nrow = nrow(fun))
    })
setMethod("E", signature(object = "AbscontCondDistribution", 
                         fun = "EuclRandMatrix", 
                         cond = "numeric"),
    function(object, fun, cond, withCond = FALSE, useApply = TRUE, ...){
        matrix(E(object, as(fun, "EuclRandVariable"), cond, withCond, 
                 useApply = useApply, ...), nrow = nrow(fun))
    })
setMethod("E", signature(object = "DiscreteCondDistribution", 
                         fun = "EuclRandMatrix", 
                         cond = "numeric"),
    function(object, fun, cond, withCond = FALSE, useApply = TRUE, ...){
        matrix(E(object, as(fun, "EuclRandVariable"), cond, withCond, 
                 useApply = useApply, ...), nrow = nrow(fun))
    })
setMethod("E", signature(object = "UnivariateCondDistribution",
                         fun = "EuclRandVarList", 
                         cond = "numeric"),
    function(object, fun, cond, withCond = FALSE, useApply = TRUE, ...){
        nrvalues <- length(fun)
        res <- vector("list", nrvalues)
        for(i in 1:nrvalues)
            res[[i]] <- E(object, fun[[i]], cond, withCond, useApply = useApply, ...)

        return(res)
    })
setMethod("E", signature(object = "AbscontCondDistribution", 
                         fun = "EuclRandVarList", 
                         cond = "numeric"),
    function(object, fun, cond, withCond = FALSE, useApply = TRUE, ...){
        nrvalues <- length(fun)
        res <- vector("list", nrvalues)
        for(i in 1:nrvalues)
            res[[i]] <- E(object, fun[[i]], cond, withCond, useApply = useApply, ...)

        return(res)
    })
setMethod("E", signature(object = "DiscreteCondDistribution", 
                         fun = "EuclRandVarList", 
                         cond = "numeric"),
    function(object, fun, cond, withCond = FALSE, useApply = TRUE, ...){
        nrvalues <- length(fun)
        res <- vector("list", nrvalues)
        for(i in 1:nrvalues)
            res[[i]] <- E(object, fun[[i]], cond, withCond, useApply = useApply, ...)

        return(res)
    })
