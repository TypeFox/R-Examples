###############################################################################
## package 'secr'
## trap.builder.R
## repeat trap layout systematically, by GRTS, or at simple random centres
## across a region
## Also make.systematic, mash(), cluster.counts(), cluster.centres()
## 2011-08-16 (full argument names)
## 2012-01-11 (cutval and signal in mash(); mxy[1:2])
## 2012-02-07 mash noise
## 2012-07-26 mash cleans out trap attributes; names sessions for list input
## 2012-09-17 mash moved to mash.r
## 2013-03-02 exclude, exclmethod arguments added
## 2013-04-20 replace deprecated overlay with over
## 2014-10-25 revamp polygon requirement for grts
## 2014-12-11 set proj4string to NA
## 2014-12-11 method = "GRTS" changed to method == "GRTS" !
###############################################################################

## spsurvey uses sp

boundarytoSPDF <- function (boundary) {
    ## build sp SpatialPolygonsDataFrame object
    ## input is 2-column matrix for a single polygon
    ## requires package sp
    Sr1 <- Polygon(boundary)
    Srs1 <- Polygons(list(Sr1), "s1")
    SpP <- SpatialPolygons(list(Srs1))
    attr <- data.frame(a = 1, row.names = "s1")
    SpatialPolygonsDataFrame(SpP, attr)
}
boundarytoSP <- function (boundary) {
    ## build sp SpatialPolygonsDataFrame object
    ## input is 2-column matrix for a single polygon
    ## requires package sp
    Sr1 <- Polygon(boundary)
    Srs1 <- Polygons(list(Sr1), "s1")
    SpatialPolygons(list(Srs1))
}
###############################################################################

trap.builder <- function (n = 10, cluster, region = NULL, frame =
    NULL, method = c("SRS", "GRTS", "all", "rank"), edgemethod =
    c("clip", "allowoverlap", "allinside"), samplefactor = 2, ranks =
    NULL, rotation = NULL, detector, exclude = NULL, exclmethod =
    c("clip", "alloutside"), plt = FALSE, add = FALSE) {

    ## region may be -
    ## matrix x,y
    ## sp SpatialPolygonsDataFrame object

    ## future: allow
    ##   minimum separation
    ##   shapefile

    #####################################################
    # 1. form polygon
    # 2. get n random origins SRS, GRTS
    # 3. translate n times; covariate for cluster of origin
    # 4. optionally clip reject edge clusters
    # 5. rbind
    #####################################################

    .local <- new.env()   ## for clusteri
    method <- match.arg(method)
    edgemethod <- match.arg(edgemethod)
    exclmethod <- match.arg(exclmethod)

    allinside <- function (xy) {
        xy <- SpatialPoints(as.matrix(xy))
        ## 2014-10-25 polygons() works with both SP and SPDF
        !any(is.na(sp::over (xy, polygons(region))))
    }

    alloutside <- function (xy) {
        xy <- SpatialPoints(as.matrix(xy))
        ## 2014-10-25 polygons() works with both SP and SPDF
        all(is.na(sp::over (xy, polygons(exclude))))
    }

    position <- function (i, cluster) {
        newtraps <- shift(cluster, origins[i,])
        if (!is.null(rotation)) {
            if (rotation<0)
                rotation <- runif(1) * 360
            newtraps <- rotate(newtraps, rotation, apply(newtraps,2,mean))
        }
        i <- .local$clusteri
        .local$clusteri <- .local$clusteri + 1    ## global update!
        clusterID(newtraps) <- factor(rep(i,nrow(newtraps)), levels=i)
        clustertrap(newtraps) <- as.numeric(polyID(newtraps))
        newtraps
    }

    ## option for single-trap clusters
    if (missing(cluster))
        cluster <- NULL
    if (is.null(cluster)) {
        if (missing(detector))
            detector <- 'multi'
        if (!(detector %in% .localstuff$pointdetectors))
            stop ("solitary detectors must be of a point detector type")
        cluster <- make.grid(nx = 1, ny = 1, detector = detector)
        edgemethod <- 'allowoverlap'
    }
    else {
        if ((attr(cluster,'detector') %in% .localstuff$polydetectors) &
            (ndetector(cluster) > 1))
            stop("clusters with multiple polygons or transects not supported")
    }

    if (method == 'all') {
        n <- nrow(frame)
    }
    if (is.null(region))
        edgemethod <- 'allowoverlap'

    SP <- inherits(region, 'SpatialPolygons')
    SPx <- inherits(exclude, 'SpatialPolygons')

    if (is.null(frame)) {
        if (is.null(region)) {
            stop ("specify at least one of 'region' or 'frame'")
        }

        if (SP) {
            ## spsurvey requires SPDF
            if ((method == 'GRTS') & (!inherits(region, 'SpatialPolygonsDataFrame'))) {
                ## 2014-12-11
                proj4string(region) <- CRS()
                attr <- data.frame(a = 1, row.names = "s1")
                region <- SpatialPolygonsDataFrame(region, attr)
            }
        }
        else {
            region <- matrix(unlist(region), ncol = 2)
            region <- rbind (region, region[1,])  # force closure of polygon
            region <- boundarytoSPDF(region)
        }
        if (SPx) {
            ## 2014-12-11
            proj4string(exclude) <- CRS()
        }
        else if (!is.null(exclude)) {
            exclude <- matrix(unlist(exclude), ncol = 2)
            exclude <- rbind (exclude, exclude[1,])  # force closure of polygon
            exclude <- boundarytoSP(exclude)
        }
        if (plt & !add) {
            plot(region)
            if (!is.null(exclude))
                plot(exclude, col='lightgrey', add=TRUE, border='lightgrey')
        }
    }
    else {
        if (plt & !add) {
            if (!is.null(region)) {
                plot(region)
            }
            else {
                MASS::eqscplot (frame, axes = F, xlab = '', ylab = '', pch = 1, cex = 0.5)
            }
            if (!is.null(exclude))
                plot(exclude, col = 'lightgrey', add = T)
        }
    }

    ntrial <- max(n * samplefactor, 5)
    ####################################
    if (method == 'SRS') {
        if (is.null(frame)) {
            origins <- coordinates(spsample(region, ntrial, type='random'))
        }
        else {
            if (ntrial > nrow(frame))
                stop ("too few rows in frame for requested sample")
            OK <- sample.int(nrow(frame), ntrial, replace = FALSE)
            origins <- as.matrix(frame[OK, ])
        }
    }
    ####################################
    else if (method == 'GRTS') {
        if (!requireNamespace ('spsurvey', quietly = TRUE))
            stop ("package 'spsurvey' required for grts in trap.builder")
        ## make a list in the format needed by grts()
        design <- list(None = list(panel = c(Panel1 = n),
            seltype = "Equal", over = ntrial))
        src <-ifelse (is.null(frame), 'sp.object', 'att.frame')
        typ <-ifelse (is.null(frame), 'area', 'finite')
        GRTS.sites <- spsurvey::grts (design = design, type.frame = typ,
            src.frame = src, sp.object = region, att.frame = frame,
            shapefile = FALSE)
        origins <- coordinates(GRTS.sites)
    }
    ####################################
    else if (method == 'all') {
        if (is.null(frame)) {
            stop ("'all' requires finite frame")
        }
        origins <- as.matrix(frame)
    }
    ####################################
    else if (method == 'rank') {
        if (is.null(frame)) {
            stop ("'rank' requires finite frame")
        }
        if (is.null(ranks)) {
            stop ("'rank' requires ranks")
        }
        nframe <- nrow(frame)
        if (nframe<n)
            stop ("not enough rows in frame for requested n")
        ranks <- ranks + runif(nframe)/(nframe+1)
        frameorder <- rev((1:nframe)[order(ranks)])
        frame <- frame[frameorder,]
        origins <- as.matrix(frame)
    }
    else
        stop ("method not recognised")

    #######################################################
    ## centre cluster on (0,0)
    if (nrow(cluster)>1)
        cxy <- apply(cluster,2,mean)
    else
        cxy <- unlist(cluster)    ## assume one detector
    cluster[,] <- sweep(cluster, MARGIN=2, FUN='-', STATS=cxy)
    #######################################################

    .local$clusteri <- 1    ## updated within position()
    if (method %in% c('all','rank')) {
        ## position all, even if we will later reject some on ranks
        traps <- lapply(1:nrow(frame), position, cluster)
    }
    else {
        traps <- lapply(1:(ntrial), position, cluster)
    }

    if (!(edgemethod %in% c('clip','allowoverlap','allinside')))
        stop ("edgemethod not recognised")

    if (!(exclmethod %in% c('clip','alloutside')))
        stop ("exclmethod not recognised")

    if ((edgemethod == 'allinside') | (exclmethod == 'alloutside')) {
        if (edgemethod %in% c('allinside')) {
            if (is.null(region))
                stop ("allinside requires 'region'")
            OK <- sapply(traps, allinside)
        }
        else OK <- length(traps)
        if (!is.null(exclude) & (exclmethod=='alloutside'))
            OK <- OK & sapply(traps, alloutside)
        if (method == 'all')
            n <- sum(OK)
        if (sum(OK) < n)
            stop ("not enough clusters inside polygon")
        traps <- traps[OK]
    }
    traps <- traps[1:n]   ## first n usable clusters

    ## convert list of clusters to flat traps object
    if (n == 1)
        traps <- traps[[1]]
    else {
        traps$renumber <- FALSE
        traps <- do.call(rbind, traps)
    }

    ## drop excluded sites, if requested
    if (edgemethod == 'clip') {
        xy <- SpatialPoints(as.matrix(traps))
        OK <- sp::over (xy, polygons(region))
        traps <- subset(traps, subset = !is.na(OK))
    }
    if (!is.null(exclude) & (exclmethod == 'clip')) {
        xy <- SpatialPoints(as.matrix(traps))
        notOK <- sp::over (xy, polygons(exclude))
        traps <- subset(traps, subset = is.na(notOK))
    }

    ## renumber clusters
    if (attr(cluster,'detector') %in% .localstuff$polydetectors) {
        npoly <- ndetector(traps)
        npercluster <- nrow(cluster)
        polyID(traps) <- factor(rep(1:npoly, rep(npercluster, npoly)))
        clustertrap(traps) <- rep(1, nrow(traps))
        clusterID(traps) <- polyID(traps)
        vertexpart <- rep(rownames(cluster), npoly)
        row.names(traps) <- paste(polyID(traps), vertexpart, sep = '.')
    }
    else {
        clusterID(traps) <- factor(as.numeric(clusterID(traps)))
        if (nrow(cluster) == 1)
            newnames <- clusterID(traps)
        else
            newnames <- paste(clusterID(traps),
                row.names(cluster)[clustertrap(traps)], sep='.')
        row.names(traps) <- newnames
    }

    ####################################
    ## optional plot
    if (plt) {
        plot(traps, add=TRUE)
        invisible(traps)
    }
    else
        traps
    ####################################
}
###############################################################################

make.systematic <- function (n, cluster, region, spacing = NULL,
    origin = NULL, ...) {

## 'cluster' is a traps object for one module
## 'region' is a rectangular survey region
## ... arguments passed to trap.builder (rotate, detector)

    SP <- inherits(region, "SpatialPolygons")
    if (SP) {
        region <- polygons(region)
        ## 2014-12-11
        proj4string(region) <- CRS()
    }
    else{
        ## convert to SpatialPolygons
        ## future: recognise & import shapefile
        region <- matrix(unlist(region), ncol = 2)
        region <- rbind (region, region[1,])  # force closure of polygon
        region <- boundarytoSP(region)
    }
    wd <- diff(bbox(region)[1,])
    ht <- diff(bbox(region)[2,])

    if (missing(cluster)) {
        ## this case is passed to trap builder for single detector placement
        ## if ... does not include detector, detector defaults to 'multi'
        cluster <- NULL
        clwd <- 0
        clht <- 0
    }
    else {
        clwd <- diff(range(cluster$x))
        clht <- diff(range(cluster$y))
    }

    wx <- clwd/2
    wy <- clht/2

    if (!is.null(spacing)) {
        rx <- spacing[1]
        ry <- ifelse(length(spacing)>1, spacing[2], rx)
        nx <- round ((wd-2*wx)/rx)
        ny <- round ((ht-2*wy)/ry)
    }
    else {
        if (length(n)>1) {
            nx <- n[1]
            ny <- n[2]
        }
        else {
            area <- sum(sapply(region@polygons, function(x) x@area))
            cell <- sqrt(area / n)
            nx <- round ((wd - 2*wx) / cell)
            ny <- round ((ht - 2*wy) / cell)
        }
        rx <- (wd - 2*wx) / nx
        ry <- (ht - 2*wy) / ny
    }
    rxy <- c(rx,ry)
    if (is.null(origin))
        origin <- runif(2) * rxy + bbox(region)[,1] + c(wx,wy)
    else {
        origin <- origin + rxy * trunc((bbox(region)[,1] - origin) / rxy)
    }
    centres <- expand.grid (
        x = seq(0, by = rx, len = nx) + origin[1],
        y = seq(0, by = ry, len = ny) + origin[2])
    centres <- SpatialPoints(as.matrix(centres))
    OK <- !is.na(sp::over (centres, region))
    centres <- coordinates(centres[OK,])
    trap.builder (cluster = cluster, frame = centres, region = region,
        method = 'all', ...)
}

###############################################################################

cluster.counts <- function (object) {
    if (!inherits(object, 'capthist'))
        stop ("requires capthist object")
    clust <- clusterID(traps(object))
    if (is.null(clust) | (length(clust) ==0) ) {
        clust <- factor(1: nrow(traps(object)))
        warning ("clusters not defined, so treating each detector as a cluster")
    }
    cl <- clust[trap(object, names = FALSE)]
    tmp <- data.frame(ID=animalID(object), cluster = cl)
    sapply(split(tmp,tmp$cluster), function(x) length(unique(x$ID)))
}
###############################################################################

cluster.centres <- function (object) {
    if (!inherits(object, 'traps'))
        stop ("requires traps object")
    clust <- clusterID(object)
    if (is.null(clust) | (length(clust) ==0) ) {
        clust <- factor(1: nrow(object))
        warning ("clusters not defined, so treating each detector as a cluster")
    }
    data.frame(x = tapply(object$x,clust,mean),
               y = tapply(object$y,clust,mean))
}
###############################################################################
