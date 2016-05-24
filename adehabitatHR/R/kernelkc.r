kernelkc <- function(tr, h, tcalc, t0, grid=40, circular=FALSE,
                     cycle=24*3600, same4all=FALSE,
                     byburst=FALSE, extent=0.5)
{
    if (!inherits(tr, "ltraj"))
        stop("tr should be of class \"ltraj\"")
    if (length(h)!=3)
        stop("h should be of length 3")
    if (!missing(t0)) {
        if (!inherits(t0, "POSIXct"))
            stop("t0 should be of class POSIXct")
    }
    if (circular) {
        if (h[3]>1)
            stop("When circular=TRUE, h should be <=1")
        h[3] <- (1-h[3])
        if (tcalc > cycle)
            stop("the time requested is beyond the cycle")
        if (missing(t0))
            stop("t0 required when circular=TRUE")
        tcalc <- unclass(tcalc)
        tcalc <- 2*pi*(tcalc%%cycle)/cycle
    } else {
        if (!inherits(tcalc, "POSIXct"))
            stop("When circular=FALSE, tcalc should be of class POSIXct")
        tcalc <- unclass(tcalc)
    }
    if (!byburst) {
        idd <- id(tr)
        tr <- lapply(unique(idd), function(i) {
            uu <- tr[id=i]
            return(do.call("rbind",uu))
        })
        names(tr) <- unique(idd)
    } else {
        names(tr) <- burst(tr)
    }
    xy <- do.call("rbind",tr)[,c("x","y")]
    xy <- xy[!is.na(xy[,1]),]
    resultats <- list()

    if (same4all) {
        if (inherits(grid, "SpatialPoints"))
            stop("when same4all is TRUE, grid should be a number")
        grid <- .makegridUD(xy, grid, extent)
    }

    for (i in 1:length(tr)) {

        if (is.list(grid)) {
            grida <- grid[names(tr)[i]][[1]]
        } else {
            grida <- grid
        }
        if (!inherits(grida, "SpatialPixels")) {
            if ((!is.numeric(grida)) | (length(grida) != 1))
                stop("grid should be a number or an object inheriting the class SpatialPixels")
            grida <- .makegridUD(xy, grida, extent)
        }
        gridded(grida) <- TRUE
        fullgrid(grida) <- TRUE
        grrw <- gridparameters(grida)
        pfs <- proj4string(grida)
        if (nrow(grrw) > 2)
            stop("grid should be defined in two dimensions")
        if (is.loaded("adehabitatMA")) {
            opteps <- adehabitatMA::adeoptions()$epsilon
        }
        else {
            opteps <- 1e-08
        }
        if ((grrw[1, 2] - grrw[2, 2]) > opteps)
            stop("the cellsize should be the same in x and y directions")
        xyg <- coordinates(grida)
        xg <- unique(xyg[, 1])
        yg <- unique(xyg[, 2])


        dft <- tr[[i]]
        dft <- dft[!is.na(dft[,1]),]
        df <- dft[, c("x", "y")]

        x <- tr[[i]]
        xy <- x[,c("x","y")]
        da <- x$date
        da <- da[!is.na(xy[,1])]
        xy <- as.matrix(xy[!is.na(xy[,1]),])

        da <- unclass(da)

        if (circular) {
            da <- (da-unclass(t0))%%cycle
            da <- (da/cycle)*2*pi
        }

        circular <- as.numeric(circular)
        xyd <- cbind(xy,da)

        toto <- .C("kernelkcr", as.double(t(xyd)), as.double(tcalc),
                   as.integer(nrow(xy)), double(length(xg)*length(yg)),
                   as.double(xg),
                   as.double(yg), as.integer(length(xg)),
                   as.integer(length(yg)),
                   as.double(h), as.integer(circular), PACKAGE="adehabitatHR")

        too <- c(t(matrix(toto[[4]], nrow=length(xg), byrow=TRUE)))
        UD <- data.frame(ud=too)
        coordinates(UD) <- expand.grid(yg,xg)[,2:1]
        gridded(UD) <- TRUE

        UD <- new("estUD", UD)
        slot(UD, "h") <- list(values=h,
                              meth="KC-specified")
        slot(UD, "vol") <- FALSE
        if (!is.na(pfs))
            proj4string(UD) <- CRS(pfs)
        resultats[[names(tr)[i]]] <- UD
    }
    names(resultats) <- names(tr)
    class(resultats) <- "estUDm"
    if (length(resultats)==1)
        resultats <- resultats[[1]]
    return(resultats)
}




exwc <- function(hv)
{
    h <- 1-hv
    x <- seq(-pi,pi,length=101)
    y <- (1 - h^2)/(2*pi*(1+h^2-2*h*cos(x)))
    plot(-50:50, y, ylab="Density", xlab="Percent of the period", ty="l")
}




kernelkcbase <- function(xyt, h, tcalc, t0, grid=40, circular=FALSE,
                         cycle=24*3600, extent=0.5)
{
    if (ncol(xyt)!=3)
        stop("xyt should have three columns")
    if (length(h)!=3)
        stop("h should be of length 3")
    if (circular) {
        if (h[3]>1)
            stop("When circular=TRUE, h should be <=1")
        h[3] <- (1-h[3])
        if (tcalc > cycle)
            stop("the time requested is beyond the cycle")
        if (missing(t0))
            stop("t0 required when circular=TRUE")
        tcalc <- unclass(tcalc)
        tcalc <- 2*pi*(tcalc%%cycle)/cycle
    } else {
        tcalc <- unclass(tcalc)
    }
    xy <- xyt[,1:2]
    gr <- grid
    resultats <- list()

    df <- xy

    if (!inherits(gr,"SpatialPixels")) {
        if (!is.numeric(gr))
            stop("grid should be an object of class SpatialPixels or a number")
        grid <- .makegridUD(xy, grid, extent)
    }
    fullgrid <- TRUE

    da <- xyt[,3]
    xy <- as.matrix(xy)

    da <- unclass(da)

    if (circular) {
        da <- (da-unclass(t0))%%cycle
        da <- (da/cycle)*2*pi
    }

    pfs <- proj4string(grid)
    xyg <- coordinates(grid)
    xg <- unique(xyg[,1])
    yg <- unique(xyg[,2])

    circular <- as.numeric(circular)
    xyd <- cbind(xy,da)

    toto <- .C("kernelkcr", as.double(t(xyd)), as.double(tcalc),
               as.integer(nrow(xy)), double(length(xg)*length(yg)),
               as.double(xg),
               as.double(yg), as.integer(length(xg)),
               as.integer(length(yg)),
               as.double(h), as.integer(circular), PACKAGE="adehabitatHR")

    too <- c(t(matrix(toto[[4]], nrow=length(xg), byrow=TRUE)))
    UD <- data.frame(ud=too)
    coordinates(UD) <- expand.grid(yg,xg)[,2:1]
    gridded(UD) <- TRUE

    UD <- new("estUD", UD)
    slot(UD, "h") <- list(values=h,
                          meth="KC-specified")
    slot(UD, "vol") <- FALSE
    if (!is.na(pfs))
        proj4string(UD) <- CRS(pfs)
    return(UD)
}


