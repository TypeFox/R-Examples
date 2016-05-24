"niche.test" <- function (x, pts, nrep = 999, o.include = TRUE, ...)
{
    ## Verifications
    if (!inherits(x, "SpatialPixelsDataFrame"))
        stop("should be an object of class SpatialPixelsDataFrame")
    gridded(x) <- TRUE
    fullgrid(x) <- FALSE
    gr <- gridparameters(x)
    if (nrow(gr) > 2)
        stop("x should be defined in two dimensions")
    if ((gr[1, 2] - gr[2, 2]) > get(".adeoptions", envir = .adehabitatMAEnv)$epsilon)
        stop("the cellsize should be the same in x and y directions")
    if (!inherits(pts, "SpatialPoints"))
        stop("should inherit from class \"SpatialPoints\"")
    pfs <- proj4string(x)

    toto <- as.matrix(join(pts, x))
    tutu <- apply(toto, 1, function(x) any(is.na(x)))
    if (sum(tutu) > 0)
        stop("points outside the study area")

    ## conversion factors -> dummy variables
    litab <- x@data
    dude <- dudi.mix(litab, scannf = FALSE)
    cw <- dude$cw
    ta <- dude$tab
    cpts <- coordinates(pts)

    ## Simulations
    res <- t(sapply(1:nrep, function(k) {
        sam <- sample(1:nrow(ta), nrow(cpts), replace=TRUE)
        tas <- ta[sam,]
        mar <- sum(apply(tas,2,mean)^2)
        tol <- sum(apply(tas,2, function(y) sum(y^2)))
        return(c(mar, tol))
    }))

    ## Observation
    coordinates(ta) <- coordinates(x)
    gridded(ta) <- TRUE
    proj4string(ta) <- CRS(pfs)

    tas <- join(pts, ta)
    mar <- sum(apply(tas,2,mean)^2)
    tol <- sum(apply(tas,2, function(y) sum(y^2)))

    ## Output
    dfxy <- as.data.frame(res)
    obs <- c(mar, tol)
    biv.test(dfxy, obs, sub = "Tests of\nmarginality\nand tolerance",
             o.include = o.include, ...)
    return(invisible(list(dfxy = dfxy, obs = obs)))
}
