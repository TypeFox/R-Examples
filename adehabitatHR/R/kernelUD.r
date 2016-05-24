## For UD storage
setClass("estUD", representation("SpatialPixelsDataFrame",
                                 h = "list", vol="logical"))


## Set as functions
as.data.frame.estUD <- function(x, row.names, optional, ...)
    as.data.frame(as(x, "SpatialPixelsDataFrame"))


setAs("estUD", "data.frame", function(from)
      as.data.frame.estUD(from))

## image of the UD
image.estUD <- function(x, ...)
    image(as(x, "SpatialGridDataFrame"), ...)

## Internal
.makegridUD <- function(xy, grid, extent)
{
    xli<-range(xy[,1])
    yli<-range(xy[,2])
    xli<-c(xli[1]-extent*abs(diff(xli)),
           xli[2]+extent*abs(diff(xli)))
    yli<-c(yli[1]-extent*abs(diff(yli)),
           yli[2]+extent*abs(diff(yli)))
    nrcol <- grid
    u <- diff(xli)
    ref <- "x"
    if (diff(yli) > diff(xli)) {
        u <- diff(yli)
        ref <- "y"
    }
    cellsize <- u/(grid - 1)
    if (ref=="x") {
        xg <- seq(xli[1], xli[2], length=nrcol)
        yg <- seq(yli[1], yli[2], by=diff(xg[1:2]))
    } else {
        yg <- seq(yli[1], yli[2], length=nrcol)
        xg <- seq(xli[1], xli[2], by=diff(yg[1:2]))
    }
    coords <- expand.grid(yg,xg)[,2:1]
    ii <- SpatialPoints(coords)
    gridded(ii) <- TRUE
    return(ii)
}


## One simple UD (UDs)
.kernelUDs <- function(xy, h="href", grid=60, hlim=c(0.1, 1.5),
                      kern = c("bivnorm", "epa"),
                      extent=0.5)
{
    ## Verifications
    kern <- match.arg(kern)
    if (!inherits(xy, "SpatialPoints"))
        stop("xy should be of class SpatialPoints")
    pfs1 <- proj4string(xy)
    xy <- coordinates(xy)
    if ((!is.numeric(h))&(h!="href")&(h!="LSCV"))
      stop("h should be numeric or equal to either \"href\" or \"LSCV\"")
    if ((h == "LSCV")&(kern == "epa"))
      stop("LSCV is not implemented with an Epanechnikov kernel")

    if (!inherits(grid, "SpatialPixels")) {
        if ((!is.numeric(grid))|(length(grid)!=1))
            stop("grid should be a number or an object inheriting the class SpatialPixels")
        grid <- .makegridUD(xy, grid, extent)
    }

    pfs2 <- proj4string(grid)
    if (!is.na(pfs2)) {
        if (!identical(pfs2,pfs1))
            stop("points and grid do not have the same proj4string")
    }

    gridded(grid) <- TRUE
    fullgrid(grid) <- TRUE
    grrw <- gridparameters(grid)
    if (nrow(grrw) > 2)
        stop("grid should be defined in two dimensions")
    if (is.loaded("adehabitatMA")) {
        opteps <- adehabitatMA::adeoptions()$epsilon
    } else {
        opteps <- 1e-08
    }
    if ((grrw[1, 2] - grrw[2, 2])> opteps)
        stop("the cellsize should be the same in x and y directions")


    xyg <- coordinates(grid)
    xg <- unique(xyg[,1])
    yg <- unique(xyg[,2])


    varx <- var(xy[, 1])
    vary <- var(xy[, 2])
    sdxy <- sqrt(0.5 * (varx + vary))
    n <- nrow(xy)
    ex <- (-1/6)
    htmp <- h
    typh <- h
    href <- sdxy * (n^ex)
    if (kern == "epa")
        href <- href * 1.77
    if (h == "href") {
        htmp <- href
    }
    if (h == "LSCV") {
        hvec <- seq(hlim[1] * href, hlim[2] * href, length = 100)
        CV <- .C("CVmise", as.integer(nrow(xy)), as.double(xy[,1]),
                 as.double(xy[, 2]), as.double(hvec), double(length(hvec)),
                 as.integer(length(hvec)), PACKAGE = "adehabitatHR")[[5]]
        htmp <- hvec[CV == min(CV)]
        if ((CV[CV == min(CV)] == CV[1]) | (CV[CV == min(CV)] ==
               CV[length(CV)]))
            warning("The algorithm did not converge \nwithin the specified range of hlim: try to increase it")
    }
    if (!is.numeric(htmp))
        stop("non convenient value for h")


    if (kern=="bivnorm") {
        toto<-.C("kernelhr", double(length(xg)*length(yg)),as.double(xg),
                 as.double(yg),
                 as.integer(length(yg)), as.integer(length(xg)),
                 as.integer(nrow(xy)), as.double(htmp),
                 as.double(xy[,1]), as.double(xy[,2]),
                 PACKAGE="adehabitatHR")
    }
    if (kern=="epa") {
        toto<-.C("kernepan", double(length(xg)*length(yg)),as.double(xg),
                 as.double(yg),
                 as.integer(length(yg)), as.integer(length(xg)),
                 as.integer(nrow(xy)), as.double(htmp),
                 as.double(xy[,1]), as.double(xy[,2]),
                 PACKAGE="adehabitatHR")
    }

    ## output
    ud <- data.frame(ud=toto[[1]])
    coordinates(ud) <- expand.grid(yg,xg)[,2:1]

    gridded(ud) <- TRUE

    if (typh=="LSCV") {
        CV<-data.frame(h=hvec, CV=CV)
        convergence<-min(CV[,2])!=CV[1,2]
        hli <-list(CV=CV, convergence=convergence, h=htmp, meth="LSCV")
    } else {
        if (typh=="href") {
            hli <- list(h=htmp, meth="href")
        } else {
            hli <- list(h=htmp, meth="specified")
        }
    }
    ud <- new("estUD", ud)
    slot(ud, "h") <- hli
    slot(ud, "vol") <- FALSE
    return(ud)
}





.boundaryk <- function(spdf, sldf, h)
{
    if (!inherits(sldf, "SpatialLines"))
        stop("sldf should be of class SpatialLinesDataFrame")
    if (!inherits(spdf, "SpatialPoints"))
        stop("spdf should be of class SpatialPoints")
    pfs1 <- proj4string(spdf)
    cool <- coordinates(sldf)
    coop <- coordinates(spdf)

    ## list of lines:
    ddff <- do.call("rbind", lapply(cool, function(x) do.call("rbind", x)))
    repp <- unlist(lapply(cool, function(x) unlist(lapply(x, nrow))))
    fac <- unlist(mapply(rep, 1:length(repp), repp))
    cool2 <- split(as.data.frame(ddff), fac)
    cool2 <- lapply(cool2, as.matrix)

    ## for each line
    ## Checks whether cos(angle) is >0
    alphaj <- unlist(lapply(cool2,
                            function(x) diff(sapply(2:nrow(x),
                                                    function(i) {
                                                        xyj <- x[i,] - x[i-1,]
                                                        atan2(xyj[2], xyj[1])
                                                    }))))
    alphaj <- sapply(alphaj,
                     function(x) ifelse(x<(-pi), x+2*pi,
                                        ifelse(x>pi, x-2*pi,x)))

    if (any(abs(alphaj)>(3.14159265359/2)))
        stop("non convenient boundary: turning angles > pi/2")

    lj <- unlist(lapply(cool2, function(xy) sapply(2:nrow(xy), function(i) {
        sqrt(sum((xy[i,]-xy[i-1,])^2))
    })))

    if (any(lj<3*h))
        stop("non convenient boundary: segments < 3*h.\n Try to decrease the smoothing parameter, or to change the boundary")

    ## for each segment, find the relocations located
    ## at a distance less that 3*h from the line
    ## distance between each reloc and each point
    res <- lapply(cool2, function(xy) {
        re <- lapply(2:nrow(xy), function(i) {
            u <- xy[i,]-xy[i-1,]
            lj <- sqrt(sum((u)^2))
            thetaj <- atan2(u[2],u[1])
            xax <- c(cos(thetaj), sin(thetaj))
            yax <- c(-sin(thetaj), cos(thetaj))
            coopc <- sweep(coop, 2, xy[i-1,])
            xp <- coopc%*%xax
            yp <- coopc%*%yax
            cons <- xp>0&xp<lj&abs(yp)<3*h
            if (sum(cons)==0)
                return(NULL)
            sig <- unique(sign(yp[cons]))
            xp <- xp[cons]
            yp <- yp[cons]
            xymirror <- cbind(xp*cos(thetaj) + yp*sin(thetaj),
                              xp*sin(thetaj) - yp*cos(thetaj))
            xymirror <- sweep(xymirror, 2, xy[i-1,], FUN="+")
            attr(xymirror, "sign") <- sig
            return(xymirror)
        })
        sig <- unique(unlist(lapply(re, function(x) attr(x, "sign"))))
        ## All the points on one side of the boundary?
        if (length(sig)>1)
            stop("Relocations beyond the boundary... is it really a boundary?")
        re <- do.call("rbind", re)
        attr(re, "sign") <- sig
        return(re)
    })
    sig <- unlist(lapply(res, function(x) ifelse(is.null(x), 1, attr(x, "sign"))))
    res <- do.call("rbind", res)


    ## for each junction, identify the concerned relocs:
    sortie <- do.call("rbind",lapply(cool2, function(xy) {
        do.call("rbind",lapply(2:(nrow(xy)-1), function(i) {
            u1 <- xy[i,] - xy[i-1,]
            u2 <- xy[i+1,] - xy[i,]
            lj <- sqrt(sum(u1^2))
            thetaj <- atan2(u1[2],u1[1])
            thetajp1 <- atan2(u2[2],u2[1])
            yax1 <- c(-sin(thetaj), cos(thetaj))
            yax2 <- c(-sin(thetajp1), cos(thetajp1))
            xax1 <- c(cos(thetaj), sin(thetaj))
            xax2 <- c(cos(thetajp1), sin(thetajp1))
            coopc <- sweep(coop, 2, xy[i,])
            coopc1 <- sweep(coop, 2, xy[i-1,])

            yp1 <- coopc1%*%yax1
            yp2 <- coopc%*%yax2
            xp1 <- coopc1%*%xax1
            xp2 <- coopc%*%xax2
            dis <- sqrt(rowSums(coopc^2))
            cons <- (xp1<lj)&(xp2>0)&(dis<3*h)
            if (sum(cons)==0)
                return(NULL)
            alphaj <- thetajp1-thetaj
            alphaj <- ifelse(alphaj<(-pi), alphaj+2*pi, alphaj)
            alphaj <- ifelse(alphaj>pi, alphaj-2*pi, alphaj)
            omegaj <- thetajp1-alphaj/2

            coopc2 <- coopc[cons,]
            xp2 <- coopc2%*%c(cos(omegaj), sin(omegaj))
            yp2 <- coopc2%*%c(-sin(omegaj), cos(omegaj))
            xymirror <- cbind(xp2*cos(omegaj)+yp2*sin(omegaj),
                              xp2*sin(omegaj)-yp2*cos(omegaj))
            xymirror <- sweep(xymirror, 2, xy[i,], FUN="+")
            return(xymirror)
        }))
    }))
    ptsadditionnels <- as.data.frame(rbind(sortie, res))
    if (nrow(ptsadditionnels)==0)
        return(NULL)
    names(ptsadditionnels) <- c("x","y")
    coordinates(ptsadditionnels) <- c("x","y")
    if (!is.na(pfs1))
        proj4string(ptsadditionnels) <- CRS(pfs1)
    attr(ptsadditionnels, "sign") <- sig
    return(ptsadditionnels)
}



## flattening the UD beyond the boundary
.fbboun <- function(estUD, sldf, sigg, h)
{
    cool <- coordinates(sldf)
    coop <- coordinates(estUD)

    ## list of lines:
    ddff <- do.call("rbind", lapply(cool, function(x) do.call("rbind", x)))
    repp <- unlist(lapply(cool, function(x) unlist(lapply(x, nrow))))
    fac <- unlist(mapply(rep, 1:length(repp), repp))
    cool2 <- split(as.data.frame(ddff), fac)
    cool2 <- lapply(cool2, as.matrix)

    if (length(sigg)!=length(cool2))
        stop("non convenient sign")

    ## identification of the points located within 6h from the boundary:
    res <- do.call("cbind",lapply(1:length(cool2), function(g) {
        xy <- cool2[[g]]
        sigg2 <- sigg[g]
        re <- do.call("cbind",lapply(2:nrow(xy), function(i) {
            u <- xy[i,]-xy[i-1,]
            lj <- sqrt(sum((u)^2))
            thetaj <- atan2(u[2],u[1])
            xax <- c(cos(thetaj), sin(thetaj))
            yax <- c(-sin(thetaj), cos(thetaj))
            coopc <- sweep(coop, 2, xy[i-1,])
            xp <- coopc%*%xax
            yp <- coopc%*%yax
            cons <- (xp>0)&(xp<lj)&(abs(yp)<6*h)&(sign(yp)==-sigg2)
            return(cons)
        }))
        rowSums(re)>0
    }))
    res <- rowSums(res)>0


    ## for each junction, identify the concerned relocs:
    ## identify the relocations within 6h from a junction
    sortie <- do.call("cbind", lapply(1:length(cool2), function(g) {
        xy <- cool2[[g]]
        sigg2 <- sigg[g]
        re <- do.call("cbind", lapply(2:(nrow(xy)-1), function(i) {
            coopc <- sweep(coop, 2, xy[i,])
            if (all(abs(coopc[,1])>6*h)|all(abs(coopc[,2])>6*h))
                return(NULL)
            dis <- sqrt(rowSums(coopc^2))
            if (all(dis>=6*h))
                return(NULL)
            u1 <- xy[i,] - xy[i-1,]
            u2 <- xy[i+1,] - xy[i,]
            lj <- sqrt(sum(u1^2))
            thetaj <- atan2(u1[2],u1[1])
            thetajp1 <- atan2(u2[2],u2[1])
            yax1 <- c(-sin(thetaj), cos(thetaj))
            yax2 <- c(-sin(thetajp1), cos(thetajp1))
            xax1 <- c(cos(thetaj), sin(thetaj))
            xax2 <- c(cos(thetajp1), sin(thetajp1))
            coopcb <- coopc[dis<6*h,]
            coopb <- coop[dis<6*h,]
            coopc1 <- sweep(coopb, 2, xy[i-1,])
            yp1 <- coopc1%*%yax1
            yp2 <- coopcb%*%yax2
            xp1 <- coopc1%*%xax1
            xp2 <- coopcb%*%xax2
            cons <- (xp1>lj)&(xp2<0)&(sign(yp1)==-sigg2)
            cons2 <- rep(FALSE, nrow(coop))
            cons2[dis<6*h] <- cons
            return(cons2)
        }))
        rowSums(re)>0
    }))

    sortie <- rowSums(sortie)>0


    df <- slot(estUD, "data")
    slot(estUD, "data") <- do.call("data.frame",lapply(df, function(x) {
        x[res|sortie] <- 0
        return(x)
    }))
    return(estUD)

}


## multiple UD estimation
kernelUD <- function (xy, h = "href", grid = 60,
                      same4all = FALSE, hlim = c(0.1, 1.5),
                      kern = c("bivnorm", "epa"), extent = 1,
                      boundary = NULL)
{
    ## Verifications
    if (!inherits(xy, "SpatialPoints"))
        stop("xy should inherit the class SpatialPoints")
    if (ncol(coordinates(xy))>2)
        stop("xy should be defined in two dimensions")
    pfs1 <- proj4string(xy)

    if (inherits(xy, "SpatialPointsDataFrame")) {
        if (ncol(xy)!=1) {
            warning("xy should contain only one column (the id of the animals)\nid ignored")
            id <- rep("a", nrow(as.data.frame(xy)))
            m <- 2
        } else {
            id <- xy[[1]]
            m <- 1
        }
    } else {
        id <- rep("a", nrow(as.data.frame(xy)))
        m <- 2
    }
    if (!is.null(boundary)) {
        if (!inherits(boundary, "SpatialLines"))
            stop("the boundary should be an object of class SpatialLines")
    }

    if (min(table(id))<5)
        stop("At least 5 relocations are required to fit an home range")
    id<-factor(id)

    xy <- as.data.frame(coordinates(xy))

    if (same4all) {
        if (inherits(grid, "SpatialPixels"))
            stop("when same4all is TRUE, grid should be a number")
        grid <- .makegridUD(xy, grid, extent)
    }

    lixy <- split(xy, id)
    res <- lapply(1:length(lixy), function(i) {
        if (is.list(grid)) {
            grida <- grid[names(lixy)[i]][[1]]
        } else {
            grida <- grid
        }
        x <- lixy[names(lixy)[i]][[1]]
        if (!is.null(boundary)) {
            bdrk <- .boundaryk(SpatialPoints(x,
                                             proj4string=CRS(as.character(pfs1))), boundary, h)
            if (!is.null(bdrk)) {
                sigg <- attr(bdrk, "sign")
                bdrk <- as.data.frame(coordinates(bdrk))
                names(bdrk) <- names(x)
                x <- rbind(x, bdrk)
            }
        }
        ud <- .kernelUDs(SpatialPoints(x, proj4string=CRS(as.character(pfs1))),
                         h = h, grid=grida, hlim=hlim,
                        kern = kern,
                        extent= extent)
        if (!is.null(boundary)){
            if (!is.null(bdrk)) {
                ud <- .fbboun(ud, boundary, sigg, h)
                slot(ud,"data")[,1] <- slot(ud,"data")[,1]/
                    (sum(slot(ud,"data")[,1])*gridparameters(ud)[1,2]^2)
            }
        }
        if (!is.na(pfs1))
            proj4string(ud) <- CRS(pfs1)
        return(ud)
    })
    names(res) <- names(lixy)
    class(res) <- "estUDm"
    if (m==2) {
        res <- res[[1]]
    }
    return(res)
}



estUDm2spixdf <- function(x)
{
    if (!inherits(x, "estUDm"))
        stop("x should be of class \"estUDm\"")
    uu <- do.call("data.frame", lapply(x, function(y) unlist(gridparameters(y))))
    ii <- all(apply(uu,1,function(y) all(y==y[1])))
    if (!ii)
        stop("this function can be used only when the same grid was used for all animals")

    res <- do.call("data.frame", lapply(x, function(y) y[[1]]))
    coordinates(res) <- coordinates(x[[1]])
    gridded(res) <- TRUE
    if (!is.na(proj4string(x[[1]])))
        proj4string(res) <- CRS(proj4string(x[[1]]))
    return(res)
}






## print method
.showestUD <- function(object)
{
    if (!inherits(object, "estUD"))
        stop("x should be an object of class estUD")
    cat("********** Utilization distribution of an Animal ************\n\n")
    if (!object@vol) {
        cat("Type: probability density\n")
    }
    cat("Smoothing parameter estimated with a ", object@h$meth,
        "parameter\n")
    cat("This object inherits from the class SpatialPixelsDataFrame.\n")
    cat("See estUD-class for more information\n\n")

    if (object@h$meth == "LSCV") {
        m <- TRUE
        if (!object@h$convergence) {
            cat("\nWARNING!! No convergence in LSCV\n")
            cat("Consider a new fit of UD using the ad hoc method for h.\n")
        }
    }

}


print.estUDm <- function(x, ...)
{
    if (!inherits(x, "estUDm"))
        stop("x should be an object of class estUDm")
    cat("********** Utilization distribution of several Animals ************\n\n")
    if (!x[[1]]@vol) {
        cat("Type: probability density\n")
    } else {
        cat("Type: volume under UD\n")
    }
    cat("Smoothing parameter estimated with a ", x[[1]]@h$meth,
        "smoothing parameter\n")
    cat("This object is a list with one component per animal.\n")
    cat("Each component is an object of class estUD\n")
    cat("See estUD-class for more information\n\n")
    if (x[[1]]@h$meth == "LSCV") {
        m <- unlist(lapply(x, function(y) y@h$convergence))
        if (any(!m)) {
            cat("\nWARNING!! No convergence in LSCV for:\n")
            print(names(m)[!m])
            cat("Consider a new fit of UD using the ad hoc method for h.\n")
        }
    }

}

## define for S4
setMethod("show", "estUD",
          function(object) .showestUD(object))


## plotLSCV:
"plotLSCV" <- function(x)
  {
      ## Verifications
      if (!inherits(x, "estUD")&!inherits(x, "estUDm"))
          stop("x should be an object of class \"estUD\" or \"estUDm\"")

      if (inherits(x, "estUD")) {
          if (x@h$meth!="LSCV")
              stop("this method can only be used when LSCV has been used for smoothing")
          plot(x@h$CV[,1], x@h$CV[,2], pch=16,
               xlab="h parameter", ylab="CV(h)", cex=0.5)
          lines(x@h$CV[,1], x@h$CV[,2])

      } else {
          if (x[[1]]@h$meth!="LSCV")
              stop("this method can only be used when LSCV has been used for smoothing")

          opar<-par(mfrow=n2mfrow(length(x)))

          ## One graph per animal
          for (i in 1:length(x)) {
              plot(x[[i]]@h$CV[,1], x[[i]]@h$CV[,2], pch=16, main=names(x)[i],
                   xlab="h parameter", ylab="CV(h)", cex=0.5)
              lines(x[[i]]@h$CV[,1], x[[i]]@h$CV[,2])
          }
          par(opar)
      }
  }


