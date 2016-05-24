"predict.enfa" <- function (object, map, nf, ...)
{
    ## Verifications
    if (!inherits(object, "enfa"))
        stop("should be an object of class \"enfa\"")
    warning("the enfa is not mathematically optimal for prediction:\n please consider the madifa instead")

    if (!inherits(map, "SpatialPixelsDataFrame"))
          stop("should be an object of class SpatialPixelsDataFrame")
      gridded(map) <- TRUE
      gr <- gridparameters(map)
      if (nrow(gr) > 2)
          stop("map should be defined in two dimensions")
      if ((gr[1, 2] - gr[2, 2])> get(".adeoptions", envir=.adehabitatMAEnv)$epsilon)
          stop("the cellsize should be the same in x and y directions")


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
    maha <- data.frame(MD=mahalanobis(Zli, center = m, cov = cov))
    coordinates(maha) <- coordinates(map)
    gridded(maha) <- TRUE

    ## Output
    return(invisible(maha))
}

