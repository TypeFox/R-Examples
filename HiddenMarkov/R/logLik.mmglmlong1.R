logLik.mmglmlong1 <- function (object, fortran=TRUE, ...){
    subjects <- unique(object$longitude)
    N <- length(subjects)
    k <- length(object$y)
    n <- k/N
    LL <- 0
    for (i in 1:N){
        a <- (object$longitude==subjects[i])
        if (object$glmfamily$family=="binomial") size <- object$size[a]
        else size <- NA
        w <- mmglm1(object$y[a], object$Pi, object$delta,
                    object$glmfamily, object$beta, object$Xdesign[a,],
                    sigma=object$sigma, nonstat=object$nonstat,
                    size=size, msg=FALSE)
        LL <- LL + logLik(w, fortran=fortran)
    }
    return(LL)
}


