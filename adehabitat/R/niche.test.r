"niche.test" <- function (kasc, points, nrep = 999, o.include = TRUE, ...)
{
    ## Verifications
    if (!inherits(kasc, "kasc"))
        stop("should be an object of class \"kasc\"")
    if (ncol(points) != 2)
        stop("points should have 2 columns")
    nrep <- nrep + 1
    toto <- join.kasc(points, kasc)
    tutu <- apply(toto, 1, function(x) any(is.na(x)))
    if (sum(tutu) > 0)
        stop("points outside the study area")


    ## conversion factors -> dummy variables
    litab <- kasc2df(kasc)
    dude <- dudi.mix(litab$tab, scannf = FALSE)
    cw <- dude$cw
    kasc <- df2kasc(dude$tab, litab$index, kasc)

    ## prepare the data for the external call
    asc <- getkasc(kasc, names(kasc)[1])
    coo <- getXYcoords(kasc)
    rc <- lapply(coo, range)
    kasc <- as.matrix(kasc)
    kasc[is.na(kasc)] <- -9999
    asc[is.na(asc)] <- -9999
    xp <- as.matrix(points)

    ## External call to the function randmargtolpts
    toto <- .C("randmargtolpts", as.double(t(xp)), as.double(rc$x),
               as.double(rc$y), as.double(t(asc)), as.double(cw),
               as.double(t(kasc)), as.double(coo$x), as.double(coo$y),
               as.double(attr(asc,"cellsize")), double(nrep),
               double(nrep), as.integer(nrep),
               as.integer(nrow(asc)), as.integer(ncol(asc)),
               as.integer(ncol(kasc)), as.integer(nrow(xp)),
               PACKAGE = "adehabitat")

    ## Output
    mar <- toto[[10]]
    tol <- toto[[11]]
    dfxy <- data.frame(marginalite = mar, tolerance = tol)[-1,]
    obs <- c(mar[1], tol[1])
    biv.test(dfxy, obs, sub = "Tests of\nmarginality\nand tolerance",
             o.include = o.include, ...)
    return(invisible(list(dfxy = dfxy, obs = obs)))
}

