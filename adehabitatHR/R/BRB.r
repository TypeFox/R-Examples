## Checks whether the boundary satisfies the required conditions
.verifboundary <- function (df, sldf, h)
{
    ## checks whether the boundary is sldf
    if (!inherits(sldf, "SpatialLines"))
        stop("sldf should be of class SpatialLinesDataFrame")

    ## Coordinates of the lines: list of the segments
    cool <- coordinates(sldf)
    coop <- as.matrix(df)
    ddff <- do.call("rbind", lapply(cool, function(x) do.call("rbind",
        x)))
    repp <- unlist(lapply(cool, function(x) unlist(lapply(x,
        nrow))))
    fac <- unlist(mapply(rep, 1:length(repp), repp))
    cool2 <- split(as.data.frame(ddff), fac)
    cool2 <- lapply(cool2, as.matrix)

    ## angles between the successive segments
    alphaj <- unlist(lapply(cool2, function(x) diff(sapply(2:nrow(x),
        function(i) {
            xyj <- x[i, ] - x[i - 1, ]
            atan2(xyj[2], xyj[1])
        }))))

    alphaj <- sapply(alphaj, function(x) ifelse(x < (-pi), x +
        2 * pi, ifelse(x > pi, x - 2 * pi, x)))
    if (any(abs(alphaj) > (3.14159265359/2)))
        stop("non convenient boundary: turning angles > pi/2")

    ## Are all the points on the same side of the line
    res <- unlist(lapply(cool2, function(xy) {
        re <- lapply(2:nrow(xy), function(i) {
            u <- xy[i, ] - xy[i - 1, ]
            lj <- sqrt(sum((u)^2))
            thetaj <- atan2(u[2], u[1])
            xax <- c(cos(thetaj), sin(thetaj))
            yax <- c(-sin(thetaj), cos(thetaj))
            coopc <- sweep(coop, 2, xy[i - 1, ])
            xp <- coopc %*% xax
            yp <- coopc %*% yax
            cons <- xp > 0 & xp < lj & abs(yp) < 3 * h
            if (sum(cons) == 0)
                return(NULL)
            sig <- unique(sign(yp[cons]))
            return(sig)
        })
        sig <- unique(unlist(re))
        if (length(sig) > 1)
            stop("Relocations beyond the boundary... is it really a boundary?")
        return(sig)
    }))
    return(res)
}


## Procedure of Benhamou and Cornelis 2010
.boundarymkde <- function (df, sldf)
{
    ## boundary: list of segments
    if (!inherits(sldf, "SpatialLines"))
        stop("sldf should be of class SpatialLinesDataFrame")
    cool <- coordinates(sldf)
    coop <- as.matrix(df[,1:2])
    ddff <- do.call("rbind", lapply(cool, function(x) do.call("rbind",
        x)))
    repp <- unlist(lapply(cool, function(x) unlist(lapply(x,
        nrow))))
    h <- max(df[,3])
    fac <- unlist(mapply(rep, 1:length(repp), repp))
    cool2 <- split(as.data.frame(ddff), fac)
    cool2 <- lapply(cool2, as.matrix)

    ## angles between successive segments
    alphaj <- unlist(lapply(cool2, function(x) diff(sapply(2:nrow(x),
        function(i) {
            xyj <- x[i, ] - x[i - 1, ]
            atan2(xyj[2], xyj[1])
        }))))
    alphaj <- sapply(alphaj, function(x) ifelse(x < (-pi), x +
        2 * pi, ifelse(x > pi, x - 2 * pi, x)))
    if (any(abs(alphaj) > (pi/2)))
        stop("non convenient boundary: turning angles > pi/2")

    ## Length of segments
    lj <- unlist(lapply(cool2, function(xy) sapply(2:nrow(xy),
        function(i) {
            sqrt(sum((xy[i, ] - xy[i - 1, ])^2))
        })))
    if (any(lj < 3 * h))
        stop("non convenient boundary: segments < 3*h.\n Try to decrease the smoothing parameter, or to change the boundary")

    ## first step: mirror of points facing the segments
    res <- lapply(cool2, function(xy) {
        re <- lapply(2:nrow(xy), function(i) {
            u <- xy[i, ] - xy[i - 1, ]
            lj <- sqrt(sum((u)^2))
            thetaj <- atan2(u[2], u[1])
            xax <- c(cos(thetaj), sin(thetaj))
            yax <- c(-sin(thetaj), cos(thetaj))
            coopc <- sweep(coop, 2, xy[i - 1, ])
            xp <- coopc %*% xax
            yp <- coopc %*% yax
            cons <- xp > 0 & xp < lj & abs(yp) < 3 * h
            if (sum(cons) == 0)
                return(NULL)
            sig <- unique(sign(yp[cons]))
            xp <- xp[cons]
            yp <- yp[cons]
            hp <- df[cons,3]
            xymirror <- cbind(xp * cos(thetaj) + yp * sin(thetaj),
                xp * sin(thetaj) - yp * cos(thetaj))
            xymirror <- sweep(xymirror, 2, xy[i - 1, ], FUN = "+")
            xymirror <- cbind(xymirror, hp)
            attr(xymirror, "sign") <- sig
            return(xymirror)
        })
        re <- do.call("rbind", re)
        return(re)
    })
    res <- do.call("rbind", res)

    ## Second step: mirror of the points located at the junction
    ## of two steps
    sortie <- do.call("rbind", lapply(cool2, function(xy) {
        do.call("rbind", lapply(2:(nrow(xy) - 1), function(i) {
            u1 <- xy[i, ] - xy[i - 1, ]
            u2 <- xy[i + 1, ] - xy[i, ]
            lj <- sqrt(sum(u1^2))
            thetaj <- atan2(u1[2], u1[1])
            thetajp1 <- atan2(u2[2], u2[1])
            yax1 <- c(-sin(thetaj), cos(thetaj))
            yax2 <- c(-sin(thetajp1), cos(thetajp1))
            xax1 <- c(cos(thetaj), sin(thetaj))
            xax2 <- c(cos(thetajp1), sin(thetajp1))
            coopc <- sweep(coop, 2, xy[i, ])
            coopc1 <- sweep(coop, 2, xy[i - 1, ])
            yp1 <- coopc1 %*% yax1
            yp2 <- coopc %*% yax2
            xp1 <- coopc1 %*% xax1
            xp2 <- coopc %*% xax2
            dis <- sqrt(rowSums(coopc^2))
            cons <- (xp1 < lj) & (xp2 > 0) & (dis < 3 * h)
            if (sum(cons) == 0)
                return(NULL)
            alphaj <- thetajp1 - thetaj
            alphaj <- ifelse(alphaj < (-pi), alphaj + 2 * pi,
                alphaj)
            alphaj <- ifelse(alphaj > pi, alphaj - 2 * pi, alphaj)
            omegaj <- thetajp1 - alphaj/2
            coopc2 <- coopc[cons, ]
            hp <- df[cons,3]
            xp2 <- coopc2 %*% c(cos(omegaj), sin(omegaj))
            yp2 <- coopc2 %*% c(-sin(omegaj), cos(omegaj))
            xymirror <- cbind(xp2 * cos(omegaj) + yp2 * sin(omegaj),
                xp2 * sin(omegaj) - yp2 * cos(omegaj))
            xymirror <- sweep(xymirror, 2, xy[i, ], FUN = "+")
            xymirror <- cbind(xymirror, hp)
            return(xymirror)
        }))
    }))
    ptsadditionnels <- as.data.frame(rbind(sortie, res))
    if (nrow(ptsadditionnels) == 0)
        return(NULL)
    names(ptsadditionnels) <- c("x", "y", "h")
    return(ptsadditionnels)
}






BRB.D <- function(ltr, Tmax=NULL, Lmin=NULL, habitat=NULL, activity=NULL)
{
    ## Checks that ltr is of class ltraj
    if (!inherits(ltr, "ltraj"))
        stop("ltr should be of class \"ltraj\"")
    if (is.null(Tmax)|is.null(Lmin))
        stop("Tmax and Lmin should be specified")
    if (Tmax < 0)
        stop("Tmax should be positive")
    if (Lmin < 0)
        stop("Lmin should be positive")

    ## If activity time is available
    if (!is.null(activity)) {
        infol <- infolocs(ltr)
        if (!(activity%in%names(infol[[1]]))) {
            stop("The activity time should be stored in the infolocs\n attribute of the object (see ?infolocs")
        }
    }

    ## If a habitat map is available
    if (!is.null(habitat)) {
        if (inherits(habitat, "SpatialGridDataFrame"))
            fullgrid(habitat) <- FALSE
        if (!inherits(habitat, "SpatialPixelsDataFrame"))
            stop("habitat should be of class \"SpatialPixelsDataFrame\"")
        if (ncol(habitat)>1)
            stop("please select only one habitat variable")

        ## get the parameters
        gp <- gridparameters(habitat)
        xll <- gp[1,1]
        yll <- gp[2,1]
        cs <- gp[1,2]
        nr <- gp[1,3]
        nc <- gp[2,3]

        ## Checks that the map covers all the relocations
        if (min(unlist(lapply(ltr, function(x) min(na.omit(x[,1])))) < xll))
            stop("relocations at the east of the eastern border of habitat map")
        if (min(unlist(lapply(ltr, function(x) min(na.omit(x[,2])))) < yll))
            stop("relocations at the south of the southern border of habitat map")
        if (max(unlist(lapply(ltr, function(x) max(na.omit(x[,1])))) > xll + cs*nr))
            stop("relocations at the west of the western border of habitat map")
        if (max(unlist(lapply(ltr, function(x) max(na.omit(x[,2])))) > yll + cs*nc))
            stop("relocations at the north of the northern border of habitat map")


        ## prepares the grid for the C function
        fullgrid(habitat) <- TRUE
        hab <- habitat[[1]]
        coo <- coordinates(habitat)
        h2 <- factor(hab[order(coo[,2], coo[,1])])
        lev <- levels(h2)
        h2 <- as.numeric(h2)
    } else {
        h2 <- 1
        xll <- 0
        yll <- 0
        cs <- 0
        nr <- 0
        nc <- 0
    }

    ## For each burst
    res <- lapply(1:length(ltr), function(i) {
        x <- ltr[[i]]
        ## keep x, y, date
        df <- x[,1:3]
        ## change date as numeric
        df[,3] <- as.numeric(df[,3])

        if (!is.null(activity)) {
            ## gets the infolocs component and the activity
            inl <- infol[[i]]
            act <- inl[,names(inl)==activity]
            if (any((na.omit(act)<0)|(na.omit(act)>1)))
                stop("activity should be a proportion.\nIt is actually >1 or <0 in these data")

            ## remove the last value (not used, often NA)
            act <- act[-length(act)]
            ## and checks whether there are missing value elsewhere
            if (any(is.na(act))) {
                stop("missing values are not allowed in the activity")
            }

            ## calculate the cumulated activity time in case of missing values
            ## for X and Y
            ## converts the proportion of activity time into
            ## an actual amount of active time
            act <- c(0, cumsum(diff(df[,3])*act))

            ## remove the missing values of df
            act <- act[!is.na(df[,1])]
            df <- df[!is.na(df[,1]),]

            ## and calculates again the proportion of activity time
            act <- diff(act)
            act <- act/diff(df[,3])

        } else {
            df <- df[!is.na(df[,1]),]
            act <- 1
        }

        ## Main calculation: Global diffusion
        tmp=.Call("CalculD", df, Tmax, Lmin, act, PACKAGE="adehabitatHR")
        DD <- data.frame(n=tmp[1], D=tmp[2])
        row.names(DD) <- c("global")
        names(DD) <- c("n","D")

        ## If habitat: D per habitat type
        if (!is.null(habitat)) {
            Dph <- .Call("calculDparhab", df, as.integer(h2-1), xll, yll, cs,
                         as.integer(nr), Lmin, as.integer(length(lev)), act, Tmax, PACKAGE="adehabitatHR")

            ## Prepare results
            Dph <- as.data.frame(Dph)
            row.names(Dph) <- lev
            names(Dph) <- c("n","D")
            DD <- rbind(DD, Dph)
            DD <- as.data.frame(DD)
        }
        return(DD)
    })
    names(res) <- burst(ltr)
    class(res) <- "DBRB"
    return(res)
}


BRB <- function(ltr, D, Tmax, Lmin, hmin, type=c("UD","ID", "RD"), radius = NULL,
                maxt = NULL, filtershort = TRUE,
                habitat = NULL, activity = NULL, grid = 200, b=FALSE,
                same4all=FALSE, extent=0.5, tau = NULL, boundary=NULL)
{
    ## Checks that the passed values are ok
    if (!inherits(ltr, "ltraj"))
        stop("ltr should be of class \"ltraj\"")
    if (is.null(Tmax)|is.null(Lmin))
        stop("Tmax and Lmin should be specified")
    if (is.null(hmin))
        stop("hmin should be specified")
    if (Tmax < 0)
        stop("Tmax should be positive")
    if (Lmin < 0)
        stop("Lmin should be positive")
    if (hmin < 0)
        stop("hmin should be positive")
    type <- match.arg(type)

    if (type!="UD") {
        if (is.null(radius))
            radius <- 3*hmin
        if (is.null(maxt))
            stop(paste("maxt should be specified when type =", type))
        if (length(radius) > 1)
            stop("Only one radius allowed in this function")
    }

    if (is.null(tau)) {
        tau <- min(unlist(lapply(ltr, function(x) na.omit(x$dt))))/10
    }
    b <- as.numeric(b)

    if (min(unlist(lapply(ltr, function(x) na.omit(x$dt))))<tau)
        stop("tau is larger than the shorter time lag between successive relocations.\nTry to decrease tau.")


    ## checks that everything is OK with D:
    if (inherits(D, "DBRB")) {

        ## correct length?
        if (length(D)!=length(ltr))
            stop("Diffusion D not of the same length as ltr")

        ## habitat passed and only one D?
        if ((nrow(D[[1]])==1)&(!is.null(habitat))) {
            warning("Only one global diffusion parameter, cannot take into account habitat")
            habitat <- NULL
        }

        ## prepares the D
        D <- lapply(D, function(x) {

            if (!is.null(habitat)) {
                ## if habitat: remove the global D and
                ## replace NaN by the global value
                zou <- x[,2]
                zou[is.na(zou)] <- zou[1]
                zou <- zou[-1]
                return(zou)
            } else {
                ## if no habitat: keep only the global D
                zou <- x[1,2]
                return(zou)
            }
        })

    } else {
        ## if no DBRB: checks that D is OK
        if (length(D)!=1)
            stop("D should be a vector of length 1 or\n an object returned by the function BRB.D")
        if (D<0)
            stop("Only positive D are allowed")

        ## prepares D for the function
        D <- rep(D, length(ltr))
        D <- lapply(1:length(ltr), function(i) D)
    }

    ## If activity is available:
    if (!is.null(activity)) {
        ## checks that the name exists in infolocs
        infol <- infolocs(ltr)
        if (!(activity%in%names(infol[[1]]))) {
            stop("The activity time should be stored in the infolocs\n attribute of the object (see ?infolocs")
        }
    }

    ## if habitat is available
    if (!is.null(habitat)) {

        ## defines the habitat map as a SpatialPixelsDataFrame
        if (inherits(habitat, "SpatialGridDataFrame"))
            fullgrid(habitat) <- FALSE
        if (!inherits(habitat, "SpatialPixelsDataFrame"))
            stop("habitat should be of class \"SpatialPixelsDataFrame\"")

        ## One habitat variable
        if (ncol(habitat)>1)
            stop("please select only one habitat variable")

        ## get the parameters of the habitat map
        gp <- gridparameters(habitat)
        xll <- gp[1,1]
        yll <- gp[2,1]
        cs <- gp[1,2]
        nr <- gp[1,3]
        nc <- gp[2,3]

        ## Checks that the map covers all the relocations
        if (min(unlist(lapply(ltr, function(x) min(na.omit(x[,1])))) < xll))
            stop("relocations at the east of the eastern border of habitat map")
        if (min(unlist(lapply(ltr, function(x) min(na.omit(x[,2])))) < yll))
            stop("relocations at the south of the southern border of habitat map")
        if (max(unlist(lapply(ltr, function(x) max(na.omit(x[,1])))) > xll + cs*nr))
            stop("relocations at the west of the western border of habitat map")
        if (max(unlist(lapply(ltr, function(x) max(na.omit(x[,2])))) > yll + cs*nc))
            stop("relocations at the north of the northern border of habitat map")


        ## prepares the grid for the C function
        fullgrid(habitat) <- TRUE
        hab <- habitat[[1]]
        coo <- coordinates(habitat)
        h2 <- factor(hab[order(coo[,2], coo[,1])])
        lev <- levels(h2)
        h2 <- as.numeric(h2)
    } else {

        ## prepares the grid for the C function
        h2 <- 1
        xll <- 0
        yll <- 0
        cs <- 0
        nr <- 0
        nc <- 0
    }

    ## Case: same grid for all animals
    if (same4all) {
        if (inherits(grid, "SpatialPoints"))
            stop("when same4all is TRUE, grid should be a number")
        xy <- do.call("rbind", ltr)[,1:2]
        xy <- xy[!is.na(xy[,1]),]
        xli <- range(xy[, 1])
        yli <- range(xy[, 2])
        xli <- c(xli[1] - extent * abs(diff(xli)), xli[2] + extent *
                 abs(diff(xli)))
        yli <- c(yli[1] - extent * abs(diff(yli)), yli[2] + extent *
                 abs(diff(yli)))
        nro <- nco <- grid
        u <- diff(xli)
        if (diff(yli) > diff(xli)) {
            u <- diff(yli)
        }
        xli <- xli[1]
        yli <- yli[1]
        csi <- u/(grid - 1)
    }

    ## case: one grid passed as arguments
    if (inherits(grid, "SpatialPixels")) {
        gp <- gridparameters(grid)
        xli <- gp[1,1]
        yli <- gp[2,1]
        nro <- gp[1,3]
        nco <- gp[2,3]
        csi <- gp[2,2]
    }


    ### Then, for each burst
    res <- lapply(1:length(ltr), function(i) {

        x <- ltr[[i]]
        ## keep x, y, date
        df <- x[,1:3]
        ## change date as numeric
        df[,3] <- as.numeric(df[,3])

        if (!is.null(activity)) {
            ## gets the infolocs component and the activity
            inl <- infol[[i]]
            act <- inl[,names(inl)==activity]
            if (any((na.omit(act)<0)|(na.omit(act)>1)))
                stop("activity should be a proportion.\nIt is actually >1 or <0 in these data")

            ## remove the last value (not used, often NA)
            act <- act[-length(act)]
            ## and checks whether there are missing value elsewhere
            if (any(is.na(act))) {
                stop("missing values are not allowed in the activity")
            }

            ## calculate the cumulated activity time in case of missing values
            ## for X and Y
            ## converts the proportion of activity time into
            ## an actual amount of active time
            act <- c(0, cumsum(diff(df[,3])*act))

            ## remove the missing values of df
            act <- act[!is.na(df[,1])]
            df <- df[!is.na(df[,1]),]

            ## and calculates again the proportion of activity time
            act <- diff(act)
            act <- act/diff(df[,3])

        } else {
            df <- df[!is.na(df[,1]),]
            act <- 1
        }

        ## prepare the grid parameters for kernel estimation
        if ((!same4all)&(!inherits(grid, "SpatialPixels"))) {
            if (is.list(grid)) {
                grida <- grid[[i]]
                if (!inherits(grida, "SpatialPixels"))
                    stop("grid should be a number or an object inheriting the class SpatialPixels")
                gp <- gridparameters(grida)
                xli <- gp[1,1]
                yli <- gp[2,1]
                nro <- gp[1,3]
                nco <- gp[2,3]
                csi <- gp[2,2]
            }  else {
                xli <- range(df[, 1])
                yli <- range(df[, 2])
                xli <- c(xli[1] - extent * abs(diff(xli)), xli[2] + extent *
                         abs(diff(xli)))
                yli <- c(yli[1] - extent * abs(diff(yli)), yli[2] + extent *
                         abs(diff(yli)))
                nro <- nco <- grid
                u <- diff(xli)
                if (diff(yli) > diff(xli)) {
                    u <- diff(yli)
                }
                xli <- xli[1]
                yli <- yli[1]
                csi <- u/(grid - 1)
            }
        }

        ## Calculates the resulting interpolation
        kk <- D[[i]]
        so <- as.data.frame(.Call("fillsegments", df, Tmax, tau,
                                  hmin, kk, Lmin, b, as.integer(h2-1),
                                  xll, yll, cs, as.integer(nr), act, as.integer(filtershort),
                                  PACKAGE="adehabitatHR"))
        names(so) <- c("x","y", "h", "t")

        ## The type of UD fillsegments should also return interpolated time
        if (type == "UD") {
            wi <- rep(1, nrow(so))
            so <- so[,1:3]
        }
        if (type == "ID") {
            wi <- 1/as.double(.Call("nvisits", so[,c(1,2,4)], radius, maxt,
                                    PACKAGE="adehabitatHR"))
            so <- so[,1:3]
        }
        if (type == "RD") {
            wi <- .Call("HRresidtime", so[,c(1,2,4)], radius, maxt, PACKAGE="adehabitatHR")
            if (all(is.na(wi)))
                warning(paste("Too large radius for burst\n",
                              "The residence time is missing for all the relocations of one burst\n"))
            so <- so[!is.na(wi),]
            wi <- 1/wi[!is.na(wi)]
            so <- so[,1:3]
        }


        ## In case of boundary: checks that it is OK and mirror the points
        if (!is.null(boundary)) {
            sigg <- .verifboundary(df[,1:2], boundary, max(so$h))
            so <- rbind(so, .boundarymkde(so, boundary))
        }

        ## Kernel smoothing
        threshh <- 3
        sor2 <- as.data.frame(.Call("mkdeb", so, xli, yli, csi, as.integer(nro),
                                    as.integer(nco), wi, as.double(threshh),
                                    PACKAGE="adehabitatHR"))
        names(sor2) <- c("x","y","dens")
        coordinates(sor2) <- c("x","y")
        gridded(sor2) <- TRUE

        ## clean the other side of the boundary
        if (!is.null(boundary)) {
            sor2 <- .fbboun(sor2, boundary, sigg, max(so$h))
        }
        sor2 <- new("estUD", sor2)
        slot(sor2, "h") <- list(values = list(hmin=hmin, D=D), interp = so,
                                meth = "BRB-specified")
        slot(sor2, "vol") <- FALSE

        return(sor2)
    })
    class(res) <- "estUDm"
    if (length(res) == 1)
        res <- res[[1]]
    return(res)
}






BRB.likD <- function(ltr, Dr=c(0.1,100),
                     Tmax = NULL, Lmin = NULL,
                     habitat = NULL, activity = NULL)
{
    if (!inherits(ltr, "ltraj"))
        stop("ltr should be of class \"ltraj\"")
    if (is.null(Tmax) | is.null(Lmin))
        stop("Tmax and Lmin should be specified")
    if (Tmax < 0)
        stop("Tmax should be positive")
    if (Lmin < 0)
        stop("Lmin should be positive")
    if (any(Dr<0))
        stop("Dr should take only positive values")
    if (length(Dr)!=2)
        stop("Dr should be of length 2")
    if (!is.null(activity)) {
        infol <- infolocs(ltr)
        if (!(activity %in% names(infol[[1]]))) {
            stop("The activity time should be stored in the infolocs\n attribute of the object (see ?infolocs")
        }
    }

    if (!is.null(habitat)) {
        fullgrid(habitat) <- FALSE
        lev <- levels(factor(habitat[[1]]))
        habitat[[1]] <- as.numeric(factor(habitat[[1]]))
        nh <- max(habitat[[1]])
        ii <- rasterize.ltraj(ltr, habitat)
    }
    if (!is.null(activity)) {
        infol <- infolocs(ltr)
    }

    val <- lapply(1:length(ltr), function(i) {
        df <- ltr[[i]]
        df[,3] <- unclass(df[,3])
        if (!is.null(activity)) {
            inl <- infol[[i]]
            act <- inl[, names(inl) == activity]
            if (any((na.omit(act) < 0) | (na.omit(act) > 1)))
                stop("activity should be a proportion.\nIt is actually >1 or <0 in these data")
            act <- act[-length(act)]
            if (any(is.na(act))) {
                stop("missing values are not allowed in the activity")
            }
            act <- c(0, cumsum(diff(df[, 3]) * act))
            act <- act[!is.na(df[, 1])]
            df <- df[!is.na(df[, 1]), ]
            act <- diff(act)
            act <- act/diff(df[, 3])
        } else {
            df <- df[!is.na(df[, 1]), ]
        }


        if (!is.null(habitat)) {
            ## get the rasterized trajectory
            tr <- ii[[i]]
            ## overlay with the map
            ov <- over(tr, geometry(habitat))
            ## get the pixels of the map
            mel <- rep(NA, length(ov))
            mel[!is.na(ov)] <- habitat[na.omit(ov),][[1]]
            ## calculate the number of unique values
            mo <- tapply(mel, tr[[1]], unique)
            ## step number
            sn <- as.numeric(names(mo))
            ## habitat type
            th <- c(-1, sapply(2:(nrow(df)-1), function(j) {
                if ((!any(sn==(j-1)))|(!any(sn==(j))))
                    return(-1)
                if ((length(mo[sn==(j-1)][[1]])>1)|(length(mo[sn==(j)][[1]])>1))
                    return(-1)
                if ((is.na(mo[sn==(j-1)][[1]]))|(is.na(mo[sn==(j)][[1]])))
                    return(-1)
                if (mo[sn==(j-1)][[1]]!=mo[sn==(j)][[1]])
                    return(-1)
                if (df[j+1,3]-df[j-1,3]>Tmax)
                    return(-1)
                d <- sqrt((df[j+1,1]-df[j-1,1])^2 + (df[j+1,2]-df[j-1,2])^2)
                l1 <- sqrt((df[j,1]-df[j-1,1])^2 + (df[j,2]-df[j-1,2])^2)
                l2 <- sqrt((df[j+1,1]-df[j,1])^2 + (df[j+1,2]-df[j,2])^2)

                if (d < Lmin)
                    return(-1)
                if (!is.null(activity)) {
                    if (act[j-1]<0.0000001)
                        return(-1)
                }
                t1 <- df[j,3]-df[j-1,3]
                t2 <- df[j+1,3]-df[j,3]
                if ((t1> 2*t2)|(t1<t2/2))
                    return(-1)
                if ((l1>2*l2)|(l1 < l2/2.0))
                    return(-1)
                return(mo[sn==(j-1)][[1]])
            }), -1)

            if (!is.null(activity)) {
                df[,3] <- c(0, cumsum(act*diff(df[,3])))
            }

            prep <- rep(0, nrow(df))
            re <- as.data.frame(do.call("rbind",lapply(1:nh, function(j) {
                prep[th==j] <- 1
                tm <- .Call("Dmv", df, Dr, prep, PACKAGE="adehabitatHR")
                if (tm[1]<4) {
                    tm[2] <- NA
                }
                return(tm)
            })))
            names(re) <- c("n","D")
            row.names(re) <- lev
        }
        prep <- c(-1, sapply(2:(nrow(df)-1), function(j) {
            if (df[j+1,3]-df[j-1,3]>Tmax)
                return(-1)
            d <- sqrt((df[j+1,1]-df[j-1,1])^2 + (df[j+1,2]-df[j-1,2])^2)
            l1 <- sqrt((df[j,1]-df[j-1,1])^2 + (df[j,2]-df[j-1,2])^2)
            l2 <- sqrt((df[j+1,1]-df[j,1])^2 + (df[j+1,2]-df[j,2])^2)
            if (d < Lmin)
                return(-1)
            if (!is.null(activity)) {
                if (act[j-1]<0.0000001)
                    return(-1)
            }
            t1 <- df[j,3]-df[j-1,3]
            t2 <- df[j+1,3]-df[j,3]
            if ((t1> 2*t2)|(t1<t2/2))
                return(-1)
            if ((l1>2*l2)|(l1 < l2/2.0))
                return(-1)
            return(1)
        }), -1)
        prep <- as.numeric(prep>0)
        tm <- .Call("Dmv", df, Dr, prep, PACKAGE="adehabitatHR")
        glob <- data.frame(n=tm[1],
                           D=tm[2])
        row.names(glob) <- c("global")
        if (!is.null(habitat))
            glob <- rbind(glob, re)
        if (any(abs(na.omit(glob[,2]-Dr[1]))<.Machine$double.eps))
            warning("The algorithm did not converge: try to decrease the lower bound of Dr")
        if (any(abs(na.omit(glob[,2]-Dr[2]))<.Machine$double.eps))
            warning("The algorithm did not converge: try to increase the upper bound of Dr")
        return(glob)
    })

    names(val) <- burst(ltr)
    class(val) <- "DBRB"
    return(val)
}
