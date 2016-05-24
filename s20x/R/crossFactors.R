crossFactors <- function (x, fac2 = NULL, ...){
    UseMethod("crossFactors")
}

crossFactors.default <- function (x, fac2 = NULL, ...)
{
    fac1 <- x
    fac1.fac2 <- paste(fac1, fac2, sep = ".")

    return(fac1.fac2)
}

crossFactors.formula <- function (formula, fac2 = NULL, data = NULL, ...){
    if(missing(formula) || length(formula) !=2){
        stop("'formula' missing or incorrect")
    }

    m <- match.call ()
    m$drop.unused.levels <- TRUE
    m[[1]] <- as.name("model.frame")
    mf <- eval(m, parent.frame())

    Terms <- attr(mf, "terms")
    fac1<-mf[,1]
    fac2<-mf[,2]

    crossFactors.default(fac1, fac2)
}
