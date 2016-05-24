"kernel.area" <- function (xy, id, h = "href", grid = 40,
                           same4all = FALSE, hlim = c(0.1, 1.5),
                           kern = "bivnorm",
                           levels = seq(20, 95, by = 5), unin = c("m", "km"),
                           unout = c("ha", "km2", "m2"), extent = 0.5)
{
    ## Computes the UD and its volume
    unin <- match.arg(unin)
    unout <- match.arg(unout)
    x <- kernelUD(xy, id, h, grid, same4all, hlim, kern, extent)
    x <- getvolumeUD(x)

    ## output objects
    area <- list()
    contours <- list()

    ## For each animal
    for (j in names(x)) {
        tmpsurf <- rep(0, length(levels))

        ## for each desired home range level
        for (i in 1:length(levels)) {
            ## identify the area corresponding to the home range of the
            ## specified level
            asc <- x[[j]]$UD
            tmp <- asc < levels[i]
            cs <- attr(asc, "cellsize")
            ## the home range size
            tmpsurf[i] <- sum(as.numeric(tmp)) * cs * cs
        }
        ## Home range size is stored for each animal
        area[[j]] <- tmpsurf
    }

    ## output
    area <- data.frame(area)
    row.names(area) <- levels
    names(area) <- names(x)
    ## output units
    if (unin == "m") {
        if (unout == "ha")
            area <- area/10000
        if (unout == "km2")
            area <- area/1e+06
    }
    if (unin == "km") {
        if (unout == "ha")
            area <- area * 100
        if (unout == "m2")
            area <- area * 1e+06
    }

    ## class "hrsize", for further plotting
    class(area) <- c("hrsize", "data.frame")
    attr(area, "units") <- unout
    return(area)
}

