###############################################################################
## package 'secr'
## make.mask.R
## 2011 10 10 transferred from methods.R
## 2012 04 10 fixed bug in ymax of bounding box
## 2012 04 11 added 'rectangular' mask type
## 2014-03-22 all polygon type with no 'traps'
###############################################################################

getCentres <- function (xy) {
    nrxy <- nrow(xy)
    if (nrxy > 1)
        (xy[-1,] + xy[-nrxy,]) / 2
    else
        xy
}

make.mask <- function (traps, buffer = 100, spacing = NULL, nx = 64, ny = 64,
    type = c("traprect", "trapbuffer", "pdot", "polygon", "clusterrect",
    "clusterbuffer", "rectangular"), poly = NULL, poly.habitat = TRUE, keep.poly = TRUE,
    check.poly = TRUE, pdotmin = 0.001, ...)
{

    type <- match.arg(type)
    if (missing(traps)) traps <- NULL
    if (ms(traps)) {         ## a list of traps objects
        if (inherits(poly, 'list') & (!is.data.frame(poly)))
            stop ("lists of polygons not implemented in 'make.mask'")
        ## 2014-09-20 now passes keep.poly and check.poly
        temp <- lapply (traps, make.mask, buffer = buffer, spacing = spacing, nx = nx, ny = ny,
                        type = type, poly = poly, poly.habitat = poly.habitat, keep.poly = TRUE,
                        check.poly = TRUE,  pdotmin = pdotmin, ...)
        class (temp) <- c('list', 'mask')
        temp
      }
    else {

        allowedType <- c('traprect','trapbuffer','polygon', 'pdot', 'clusterrect',
                         'clusterbuffer', 'rectangular')
        if (! (type %in% allowedType))
            stop ("mask type must be one of ",
                  paste(sapply(allowedType, dQuote), collapse=","))
        dots <- match.call(expand.dots = FALSE)$...
        if ((length(dots)==0) & (type == 'pdot'))
            warning ("no detection parameters supplied; using defaults")

        buff <- c(-buffer,+buffer)

        if (!is.null(poly)) {
            SPDF <- inherits(poly, "SpatialPolygonsDataFrame")
            if (!SPDF) {
                poly <- matrix(unlist(poly), ncol = 2)
                poly <- rbind (poly, poly[1,])  # force closure of poly
            }
        }

        ## if (is.null(traps))
        ##    type <- 'rectangular'
        ## replaced 2014-03-22
        if (is.null(traps)) check.poly <- FALSE
        if (is.null(traps) & !(type == 'polygon'))
            type <- 'rectangular'
        if (type == 'rectangular') {
            if (is.null(spacing))
                stop ("require spacing for rectangular mask")
            xl <- c(0, spacing * nx)
            yl <- c(0, spacing * ny)
        }
        else if (type=='polygon') {
            if (is.null(poly))
                stop ("mask polygon must be supplied")
            if (!poly.habitat)
                stop ("type = 'polygon' not compatible with nonhabitat")
            if (SPDF) {
                xl <- poly@bbox[1,]
                yl <- poly@bbox[2,]
            }
            else {
                xl <- range(poly[,1])
                yl <- range(poly[,2])
            }
        }
        else {
            xl <- range(traps$x) + buff
            yl <- range(traps$y) + buff
        }

        if (is.null(spacing)) spacing <- diff(xl) / nx

        if (type %in% c('clusterrect', 'clusterbuffer')) {
            ID <- clusterID(traps)
            meanx <- unique(tapply(traps$x, ID, mean))
            meany <- unique(tapply(traps$y, ID, mean))
            cluster <- subset(traps, subset = clusterID(traps)==1) ## extract a single cluster
            ## assume identical wx, wy are half-width and half-height of a box
            ## including the cluster and the rectangular buffer
            wx <- diff(range(cluster$x)) / 2 + buffer
            wy <- diff(range(cluster$y)) / 2 + buffer
            wx <- round(wx/spacing) * spacing   ## to make symmetrical
            wy <- round(wy/spacing) * spacing   ## to make symmetrical
            dx <- seq(-wx,wx,spacing)
            dy <- seq(-wy,wy,spacing)
            x <- as.numeric(outer(FUN='+', dx, meanx))
            y <- as.numeric(outer(FUN='+', dy, meany))
        }
        else {
            x <- seq(xl[1] + spacing/2, xl[2], spacing)
            y <- seq(yl[1] + spacing/2, yl[2], spacing)

        }

        mask   <- expand.grid (x=x, y=y)
        attr(mask,'out.attrs') <- NULL   ## added 2009 07 03

        if (type=='trapbuffer') {
            ## appropriate convex buffer 2011-01-22
            ## (this re-use of nx may not be appropriate)
            if (detector(traps) %in% c('polygon','polygonX')) {
               temp <- buffer.contour(traps, buffer = buffer, nx = nx,
                                      convex = T, plt = F)
               OK <- array(dim=c(length(x), length(y), length(temp)))
               for (i in 1:length(temp))
                  OK[,,i] <- pointsInPolygon(mask, temp[[i]][,c('x','y')])
               OK <- apply(OK, 1:2, any)
               mask <- mask[OK,,drop=F]
           }
            else
                mask <- mask[distancetotrap(mask, traps) <= buffer,]
        }

        if (type=='clusterbuffer') {
            mask <- mask[distancetotrap(mask, traps) <= buffer,]
        }

        if (type=='pdot') {
            OK <- pdot(mask, traps = traps, ...) > pdotmin
            edge <- function (a,b) any (abs(a-b) < (spacing))
            mask <- mask[OK,]
            attr(mask,'pdotmin') <- pdotmin   # save nominal threshold
            if (edge(mask[,1],xl[1]) |
                edge(mask[,1],xl[2]) |
                edge(mask[,2],yl[1]) |
                edge(mask[,2],yl[2]))
            warning ("'pdot' mask may have been truncated; ",
                     "possibly increase buffer")
        }
        if (!is.null(poly)) {
            if (poly.habitat) {
                inpoly <- pointsInPolygon(mask, poly)
                mask <- mask[inpoly,]
                if (check.poly) {
                    if (any (!pointsInPolygon(traps, poly)))
                    warning ("some traps are not inside habitat polygon")
                }
            }
            else {
                mask <- mask[!pointsInPolygon(mask, poly),]
                if (check.poly)
                    if (any (pointsInPolygon(traps, poly)))
                    warning ("some traps are inside non-habitat polygon")
            }
            if (keep.poly) {
                attr(mask, 'polygon') <- poly   # save
                attr(mask, 'poly.habitat') <- poly.habitat   # save
            }
        }

        ## 2012 04 10 revised limits for bounding box
        xl <- range(mask$x) + spacing/2 * c(-1,1)
        yl <- range(mask$y) + spacing/2 * c(-1,1)

        attr(mask,'type')        <- type
        attr(mask,'meanSD')      <- getMeanSD (mask)
        attr(mask,'area')        <- spacing^2 * 0.0001
        attr(mask,'spacing')     <- spacing
        attr(mask,'boundingbox') <- expand.grid(x=xl,y=yl)[c(1,2,4,3),]
        class(mask)  <- c('mask', 'data.frame')

        mask
    }
}
###############################################################################
