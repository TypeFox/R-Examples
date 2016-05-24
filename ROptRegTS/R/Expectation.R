.condE <- function(x, D1, fct, ...){ 
    FUN <- function(x, D1, fct, ...){ 
        return(E(D1, fct, cond = x, withCond = TRUE, ...))
    }
    return(apply(matrix(x, ncol = dimension(Range(cond(D1)))), 1, FUN, D1 = D1, fct = fct, ...))
}

setMethod("E", signature(object = "L2RegTypeFamily", 
                         fun = "EuclRandVariable", 
                         cond = "missing"),
    function(object, fun){
        fct <- function(x, cond, f1){ f1(cbind(t(cond),x)) }
        
        res <- numeric(length(fun))
        for(i in 1:length(fun)){
            res[i] <- E(object@RegDistr, .condE, D1 = object@distribution, 
                        fct = fct, f1 = fun@Map[[i]])
        }
        
        return(res)
    })
setMethod("E", signature(object = "L2RegTypeFamily", 
                         fun = "EuclRandMatrix", 
                         cond = "missing"),
    function(object, fun){
        matrix(E(object, as(fun, "EuclRandVariable")), nrow = nrow(fun))
    })
setMethod("E", signature(object = "L2RegTypeFamily", 
                         fun = "EuclRandVarList", 
                         cond = "missing"),
    function(object, fun){
        nrvalues <- length(fun)
        res <- vector("list", nrvalues)
        for(i in 1:nrvalues) res[[i]] <- E(object, fun[[i]])

        return(res)
    })
