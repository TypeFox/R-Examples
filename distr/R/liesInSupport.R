setMethod("liesInSupport", signature(object = "DiscreteDistribution",
                                     x = "numeric"),
    function(object, x){ x %in% support(object) })

setMethod("liesInSupport", signature(object = "DiscreteDistribution",
                                     x = "matrix"),
    function(object, x){ 
        if(ncol(x) != 1)
            stop("'x' has wrong dimension")
        else
            as.vector(x) %in% support(object) 
    })

setMethod("liesInSupport", signature(object = "AbscontDistribution",
                                     x = "numeric"),
    function(object, x){ 
        if(!is.nan(q(object)(0)))
            low <- q(object)(0)
        else
            low <- q(object)(10*.Machine$double.eps)
        if(!is.nan(q(object)(1)))
            upp <- q(object)(1)
        else
            upp <- q(object)(1-10*.Machine$double.eps)

        (x >= low)&(x <= upp)
    })

setMethod("liesInSupport", signature(object = "AbscontDistribution",
                                     x = "matrix"),
    function(object, x){ 
        if(ncol(x) != 1)
            stop("'x' has wrong dimension")
        else
            liesInSupport(as.vector(x), object)
    })

