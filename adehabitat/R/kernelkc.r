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
    gr <- grid
    resultats <- list()

    if (is.list(gr)) {
        if (is.null(names(gr)))
            stop("when grid is a list, it should have named elements")
        nn <- names(gr)
        lev <- names(tr)
        if (length(lev) != length(nn))
            stop("the length of the grid list should be equal to the number of levels of id/burst")
        if (!all(lev %in% nn))
            stop("some levels of id/bursts do not have corresponding grids")
    }
    if (same4all) {
        if (!is.list(gr)) {

            if (length(as.vector(gr)) == 1) {
                if (!is.numeric(gr))
                    stop("grid should be an object of class asc or a number")
                xli <- range(xy[, 1])
                yli <- range(xy[, 2])
                xli <- c(xli[1] - extent * abs(xli[2] - xli[1]),
                         xli[2] + extent * abs(xli[2] - xli[1]))
                yli <- c(yli[1] - extent * abs(yli[2] - yli[1]),
                         yli[2] + extent * abs(yli[2] - yli[1]))
                xygg <- data.frame(x = xli, y = yli)
                grid <- ascgen(xygg, nrcol = grid)
                cellsize <- attr(grid, "cellsize")
                lx <- nrow(grid) * cellsize
                ly <- ncol(grid) * cellsize
                ref <- lx
                if (ly > lx)
                    ref <- ly
                xll <- attr(grid, "xll")
                yll <- attr(grid, "yll")
                xll <- xll - lx * extent
                yll <- yll - ly * extent
                arajlig <- ceiling((lx * extent)/cellsize)
                arajcol <- ceiling((ly * extent)/cellsize)
                mrajlig <- matrix(0, ncol = ncol(grid), nrow = arajlig)
                grid <- rbind(mrajlig, grid, mrajlig)
                mrajcol <- matrix(0, ncol = arajcol, nrow = nrow(grid))
                grid <- cbind(mrajcol, grid, mrajcol)
                attr(grid, "xll") <- xll
                attr(grid, "yll") <- yll
                attr(grid, "cellsize") <- cellsize
                attr(grid, "type") <- "numeric"
                class(grid) <- "asc"
            }
        }  else {
            stop("when same4all=TRUE, a list of grid cannot be passed as \"grid\"")
        }
    }

    for (i in 1:length(tr)) {

        dft <- tr[[i]]
        dft <- dft[!is.na(dft[,1]),]
        df <- dft[, c("x", "y")]
        if (!is.list(gr)) {
            if (length(as.vector(gr)) == 1) {
                if (!is.numeric(gr))
                    stop("grid should be an object of class asc or a number")
                if (!same4all) {
                    grid <- matrix(0, ncol = gr, nrow = gr)
                    rgx <- range(df[, 1])
                    rgy <- range(df[, 2])
                    lx <- rgx[2] - rgx[1]
                    ly <- rgy[2] - rgy[1]
                    ref <- lx
                    if (ly > lx)
                        ref <- ly
                    xll <- rgx[1]
                    yll <- rgy[1]
                    cellsize <- ref/ncol(grid)
                    xll <- xll - extent * lx
                    yll <- yll - extent * ly
                    arajlig <- ceiling((lx * extent)/cellsize)
                    arajcol <- ceiling((ly * extent)/cellsize)
                    mrajlig <- matrix(0, ncol = ncol(grid), nrow = arajlig)
                    grid <- rbind(mrajlig, grid, mrajlig)
                    mrajcol <- matrix(0, ncol = arajcol, nrow = nrow(grid))
                    grid <- cbind(mrajcol, grid, mrajcol)
                    attr(grid, "xll") <- xll
                    attr(grid, "yll") <- yll
                    attr(grid, "cellsize") <- cellsize
                    attr(grid, "type") <- "numeric"
                    class(grid) <- "asc"
                }
            }
        } else {
            grid <- gr[[names(lixy)[i]]]
        }

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


        lixy <- getXYcoords(grid)
        grrr <- grid
        grid[is.na(grid)] <- 0

        xg <- lixy$x
        yg <- lixy$y
        class(grid) <- "matrix"

        circular <- as.numeric(circular)
        xyd <- cbind(xy,da)

        toto <- .C("kernelkcr", as.double(t(xyd)), as.double(tcalc),
                   as.integer(nrow(xy)), as.double(t(grid)), as.double(xg),
                   as.double(yg), as.integer(nrow(grid)),
                   as.integer(ncol(grid)),
                   as.double(h), as.integer(circular), PACKAGE="adehabitat")

        UD <- matrix(toto[[4]], nrow = nrow(grrr), byrow = TRUE)
        UD <- UD/(sum(UD)*(attr(grrr, "cellsize")^2))
        UD <- getascattr(grrr, UD)
        attr(UD,"UD") <- "simple"
        resultats[[i]] <- UD
        grid <- grrr
    }
    names(resultats) <- names(tr)
    class(resultats) <- "kernelkco"
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

    if (!inherits(gr,"asc")) {
        if (!is.numeric(gr))
            stop("grid should be an object of class asc or a number")

        grid <- matrix(0, ncol = gr, nrow = gr)
        rgx <- range(df[, 1])
        rgy <- range(df[, 2])
        lx <- rgx[2] - rgx[1]
        ly <- rgy[2] - rgy[1]
        ref <- lx
        if (ly > lx)
            ref <- ly
        xll <- rgx[1]
        yll <- rgy[1]
        cellsize <- ref/ncol(grid)
        xll <- xll - extent * lx
        yll <- yll - extent * ly
        arajlig <- ceiling((lx * extent)/cellsize)
        arajcol <- ceiling((ly * extent)/cellsize)
        mrajlig <- matrix(0, ncol = ncol(grid), nrow = arajlig)
        grid <- rbind(mrajlig, grid, mrajlig)
        mrajcol <- matrix(0, ncol = arajcol, nrow = nrow(grid))
        grid <- cbind(mrajcol, grid, mrajcol)
        attr(grid, "xll") <- xll
        attr(grid, "yll") <- yll
        attr(grid, "cellsize") <- cellsize
        attr(grid, "type") <- "numeric"
        class(grid) <- "asc"
    }

    da <- xyt[,3]
    xy <- as.matrix(xy)

    da <- unclass(da)

    if (circular) {
        da <- (da-unclass(t0))%%cycle
        da <- (da/cycle)*2*pi
    }

    lixy <- getXYcoords(grid)
    grrr <- grid
    grid[is.na(grid)] <- 0

    xg <- lixy$x
    yg <- lixy$y
    class(grid) <- "matrix"

    circular <- as.numeric(circular)
    xyd <- cbind(xy,da)

    toto <- .C("kernelkcr", as.double(t(xyd)), as.double(tcalc),
               as.integer(nrow(xy)), as.double(t(grid)), as.double(xg),
               as.double(yg), as.integer(nrow(grid)),
               as.integer(ncol(grid)),
               as.double(h), as.integer(circular), PACKAGE="adehabitat")


    UD <- matrix(toto[[4]], nrow = nrow(grrr), byrow = TRUE)
    UD <- UD/(sum(UD)*(attr(grrr, "cellsize")^2))
    UD <- getascattr(grrr, UD)
    attr(UD, "UD") <- "simple"
    return(UD)
}



getvolumeUDs <- function(asc)
{
    if (!inherits(asc,"asc"))
        stop("asc should be an UD of class \"asc\"")
    if (is.null(attr(asc,"UD")))
        stop("only UD can be passed to this function")
    if (attr(asc,"UD")!="simple")
        stop("only UD can be passed to this function")

    cs <- attr(asc, "cellsize")
    asc <- asc/(sum(asc)*cs*cs)
    v <- .C("calcvolume", as.double(t(asc)),
            as.integer(ncol(asc)),
            as.integer(nrow(asc)),
            as.double(cs),
            PACKAGE = "adehabitat")[[1]]
    index <- 1:length(v)
    vord <- v[order(v, decreasing = TRUE)]
    indord <- index[order(v, decreasing = TRUE)]
    vsu <- cumsum(vord)
    vreord <- vsu[order(indord)] * 100
    u <- matrix(vreord, ncol = ncol(asc), byrow = TRUE)
    asc <- getascattr(asc, u)
    attr(asc, "UD") <- "volume"
    return(asc)
}


getvolumeUDk <- function(kc)
{
    if (!inherits(kc,"kernelkco"))
        stop("kc should be an object of class \"kernelkco\"")
    uu <- list()
    ret <- lapply(kc, function(asc) {
        getvolumeUDs(asc)
    })
    class(ret) <- "kernelkcv"
    return(ret)
}


getverticeshrk <- function(kc, lev=95)
{
    if (inherits(kc, "kernelkco"))
        kc <- getvolumeUDk(kc)
    if (!inherits(kc, "kernelkcv"))
        stop("kc should be of class kernelkcv or kernelkco")
    if (length(lev)>1)
        stop("lev should be of length 1")
      contour<-list()

      ## for each animal
    for (i in 1:length(kc)) {

        ## gets the UD and keep areas upper than lev
        ud<-kc[[i]]

        ## gets the contour of the connected features
        xyl <- getXYcoords(ud)
        re <- contourLines(x = xyl$x,
                           y = xyl$y,
                           ud, nlevels = 1,
                           levels = lev)
        so <- do.call("rbind", lapply(1:length(re), function(i) {
            so <- data.frame(fac=rep(i, length(re[[i]]$x)),
                             x=re[[i]]$x,
                             y=re[[i]]$y)
            return(so)
        }))
        so[,1] <- factor(so[,1])

        contour[[i]]<-as.area(so)
    }
    ## output of class "kver"
    names(contour)<-names(kc)
    class(contour) <- "kver"
    return(contour)
}




getverticeshrs <- function(ascv, lev=95)
{
    if (!inherits(ascv, "asc"))
        stop("asc should be of class asc")
    if (length(lev)>1)
        stop("lev should be of length 1")
    if (is.null(attr(ascv,"UD")))
        stop("only UD can be passed to this function")
    if (attr(ascv,"UD")!="volume")
        ascv <- getvolumeUDs(ascv)

    ## gets the contour of the connected features
    xyl <- getXYcoords(ascv)
    re <- contourLines(x = xyl$x,
                       y = xyl$y,
                       ascv, nlevels = 1,
                       levels = lev)
    so <- do.call("rbind", lapply(1:length(re), function(i) {
        so <- data.frame(fac=rep(i, length(re[[i]]$x)),
                         x=re[[i]]$x,
                         y=re[[i]]$y)
        return(so)
    }))
    so[,1] <- factor(so[,1])

    so <- as.area(so)
    return(so)
}
