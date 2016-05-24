
### NA handling for objects of class `ModelEnv'

complete.cases.ModelEnv <- function(x) {
    
    do.call("complete.cases", as.data.frame(lapply(ls(x@env), function(o) x@get(o))))

}

# setGeneric("na.fail", useAsDefault = na.fail)

setMethod("na.fail", signature = "ModelEnv", definition = function(object, ...) {

    cc <- complete.cases.ModelEnv(object)
    if (!all(cc)) return(FALSE)
    return(object)
})

# setGeneric("na.pass", useAsDefault = na.pass)

setMethod("na.pass", signature = "ModelEnv", definition = function(object, ...) {  

    return(object)

})

# setGeneric("na.omit", useAsDefault = na.omit)

setMethod("na.omit", signature = "ModelEnv", definition = function(object, ...) {

    cc <- complete.cases.ModelEnv(object)
    if (!all(cc)) return(subset(object, cc, ...))
    return(object)

})

