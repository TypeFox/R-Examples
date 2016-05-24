"predict.enfa" <- function (object, index, attr, nf, ...)
{
    ## Verifications
    if (!inherits(object, "enfa"))
        stop("should be an object of class \"enfa\"")
    warning("the enfa is not mathematically optimal for prediction:\n please consider the madifa instead")

    ## The number of axes of specialization for the prediction
    if ((missing(nf)) || (nf > object$nf))
        nf <- object$nf

    ## ... and also keeps the marginality axis
    Zli <- object$li[, 1:(nf + 1)]

    ## The Mahalanobis distances computed on these axes
    f1 <- function(x) rep(x, object$pr)
    Sli <- apply(Zli, 2, f1)
    m <- apply(Sli, 2, mean)
    cov <- t(as.matrix(Sli)) %*% as.matrix(Sli)/nrow(Sli)
    maha <- mahalanobis(Zli, center = m, cov = cov)

    ## Output
    map <- getkasc(df2kasc(data.frame(toto = maha, tutu = maha),
                           index, attr), "toto")
    return(invisible(map))
}

