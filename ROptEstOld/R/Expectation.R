setMethod("E", signature(object = "L2ParamFamily", 
                         fun = "EuclRandVariable", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE){
        return(E(object = object@distribution, fun = fun, useApply = useApply))
    })
setMethod("E", signature(object = "L2ParamFamily", 
                         fun = "EuclRandMatrix", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE){
        matrix(E(object = object, fun = as(fun, "EuclRandVariable"), 
                 useApply = useApply), nrow = nrow(fun))
    })
setMethod("E", signature(object = "L2ParamFamily", 
                         fun = "EuclRandVarList", 
                         cond = "missing"),
    function(object, fun, useApply = TRUE){
        nrvalues <- length(fun)
        res <- vector("list", nrvalues)
        for(i in 1:nrvalues) res[[i]] <- E(object = object, fun = fun[[i]], 
                                           useApply = useApply)

        return(res)
    })
setMethod("E", signature(object = "Gumbel",
                         fun = "missing",
                         cond = "missing"),
    function(object, low = NULL, upp = NULL, ...){a <- loc(object); b <- scale(object)
    if(is.null(low) && is.null(upp))
           return(a- EULERMASCHERONICONSTANT * b)
    else
        return(E(as(object,"AbscontDistribution"), low=low, upp=upp, ...))
    })
## http://mathworld.wolfram.com/GumbelDistribution.html
