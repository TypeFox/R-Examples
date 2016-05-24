## "Math" group
setMethod("Math", signature(x = "EuclRandVariable"),
    function(x){
        nrvalues <- length(x)
        map <- vector("list", nrvalues)

        fct1 <- NULL; f <- function(x){};
        for(i in 1:nrvalues){
            map[[i]] <- function(x){ f1 <- fct1; f(f1(x)) }
            body(map[[i]]) <- substitute({ f1 <- fct1; f(f1(x)) },
                                    list(f = as.name(.Generic), fct1 = x@Map[[i]]))
        }
        
        x@Map <- map
        return(x)
    })
setMethod("Math", signature(x = "EuclRandMatrix"),
    function(x){
        x@Map <- Map(callGeneric(as(x, "EuclRandVariable")))

        return(x)
    })
setMethod("Math", signature(x = "EuclRandVarList"),
    function(x){
        for(i in 1:length(x))
            x[[i]] <- callGeneric(x[[i]])

        return(x)
    })
