
raster <- function(theta, pred, fam, root, famlist = fam.default())
{
    stopifnot(is.matrix(theta))
    stopifnot(is.numeric(theta))
    nind <- nrow(theta)
    nnode <- ncol(theta)
    storage.mode(theta) <- "double"

    stopifnot(is.numeric(pred))
    stopifnot(length(pred) == nnode)
    stopifnot(all(pred == as.integer(pred)))
    stopifnot(all(pred < seq(along = pred)))

    stopifnot(is.numeric(fam))
    stopifnot(length(fam) == nnode)
    stopifnot(all(fam == as.integer(fam)))
    stopifnot(is.element(fam, seq(along = famlist)))

    stopifnot(is.matrix(root))
    stopifnot(is.numeric(root))
    stopifnot(nind == nrow(root))
    stopifnot(nnode == ncol(root))
    storage.mode(root) <- "double"

    setfam(famlist)

    result <- .C("aster_simulate_data",
        nind = as.integer(nind),
        nnode = as.integer(nnode),
        pred = as.integer(pred),
        fam = as.integer(fam),
        theta = theta,
        root = root,
        x = matrix(as.double(0), nind, nnode),
        PACKAGE = "aster")$x

    clearfam()
    return(result)
}

