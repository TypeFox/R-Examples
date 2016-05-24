mex <- function(data, which, mth, mqu, dqu, margins="laplace",constrain=TRUE,v=10,
                penalty="gaussian", maxit=10000,
                trace=0, verbose=FALSE, priorParameters=NULL){

    theCall <- match.call()

    if(is.null(colnames(data))){
        colnames(data) <- paste(rep("Column",ncol(data)),1:ncol(data),sep="")
    }

    if (missing(which)){
        which <- colnames(data)[1]
        cat("which not given. Conditioning on", which, "\n")
    }

    if (missing(mth)){
        if(length(mqu) == 1){
            mqu <- rep(mqu, ncol(data))
    }
    if( length(mqu) != ncol(data)){
      stop("mqu must be length 1 or length equal to the dimension of the data")
    }
    mth <- unlist(lapply(1:length(mqu), function(i, data, p){
                                          quantile(data[, i], prob=p[i])
                                        }, p=mqu, data=data))
    }

    res1 <- migpd(data=data, mth=mth, penalty=penalty,
                  maxit=maxit, trace=trace, verbose=verbose,
                  priorParameters=priorParameters)

    res2 <- mexDependence(x= res1, which=which, dqu=dqu, margins=margins, constrain=constrain, v=v)
    res2$call <- theCall
    res2
}

print.mex <- function(x, ...){
    print(x$call, ...)
    cat("\n\nMarginal models:\n")
    summary(x[[1]])
    cat("\nDependence model:\n\n")
    print(x[[2]])
    invisible()
}

summary.mex <- function(object, ...){
    print(object, ...)
    invisible(coef(object))
}

coefficients.mex <- function(object, ...){
    res1 <- coef(object[[1]])
    res2 <- coef(object[[2]]) # uses native coef method
    list(margins=res1, dependence=res2)
}

coef.mex <- coefficients.mex
