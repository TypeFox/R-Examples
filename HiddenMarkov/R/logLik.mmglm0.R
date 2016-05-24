"logLik.mmglm0" <- function(object, fortran=TRUE, ...){
    object <- as.dthmm(object)
    return(logLik.dthmm(object, fortran=fortran))
}

