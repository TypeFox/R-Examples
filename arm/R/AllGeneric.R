

#setGeneric("coef")
#setGeneric("print")
#setGeneric("fitted")

#setGeneric("extractAIC")

if (!isGeneric("coefplot")) {
    setGeneric("coefplot",
               function(object, ...)
               standardGeneric("coefplot"))
}


if (!isGeneric("display")) {
    setGeneric("display",
               function(object, ...)
               standardGeneric("display"))
}


if (!isGeneric("sim")) {
    setGeneric("sim",
               function(object, ...)
               standardGeneric("sim"))
}

sigma.hat <- function(object,...){
    UseMethod("sigma.hat")
}


if (!isGeneric("se.coef")) {
    setGeneric("se.coef",
               function(object, ...)
               standardGeneric("se.coef"))
}


if (!isGeneric("mcsamp")) {
    setGeneric("mcsamp",
               function(object, ...)
               standardGeneric("mcsamp"))
}



if (!isGeneric("standardize")) {
    setGeneric("standardize",
               function(object, ...)
               standardGeneric("standardize"))
}



#if (!isGeneric("terms.bayes")) {
#    setGeneric("terms.bayes",
#               function(x, ...)
#               standardGeneric("terms.bayes"))
#}



if (!isGeneric("traceplot")) {
    setGeneric("traceplot",
               function(x, ...)
               standardGeneric("traceplot"),
               useAsDefault = function(x, ...) coda::traceplot(x, ...))
}
