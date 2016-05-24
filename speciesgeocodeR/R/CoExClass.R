CoExClass <- function(x) {
    if (class(x) == "spgeoOUT") {
        pp <- .CoExClassH(x$spec_table)
        x$coexistence_classified <- pp
        return(x)
    } else {
        stop("function is only defined for class \"SpgeoOUT\". \n  See .CoExClassH() for single \"data.frames\"")
    }
} 
