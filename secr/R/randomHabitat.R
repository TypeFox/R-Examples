## Simulated categorical landscape
## Modified Random Cluster method of Saura and Martinez-Millan 2000
## Landscape Ecology 15: 661-678

## Murray Efford 2012-04-09,10,11,12;
## 2013-07-09 replaced call to obsolete raster function 'adjacency'
## 2014-08-25 'raster' now in Imports, so do not 'require'
## Requires 'raster' and 'igraph0' packages:

## Robert J. Hijmans & Jacob van Etten (2011). raster: Geographic
##  analysis and modeling with raster data. R package version 1.9-33.
##  http://CRAN.R-project.org/package=raster

## Csardi G, Nepusz T: The igraph software package for complex network
##  research, InterJournal, Complex Systems 1695. 2006.
##  http://igraph.sf.net

## Restrictions
##   adjacency for step 'B' is based on a 1-step rook move (direction = 4)
##   (this yields an isotropic result)

randomHabitat <- function (mask, p = 0.5, A = 0.5, directions = 4, minpatch = 1,
                           drop = TRUE, covname = 'habitat', plt = FALSE) {

    if (ms(mask)) {
        ## allow for possibility that 'mask' is a list of masks
        temp <- lapply(mask, randomHabitat, p = p, A = A, directions = directions,
                       minpatch = minpatch, drop = drop, covname = covname, plt = plt)
        class(temp) <- c('list','mask')
        temp
    }
    else {

        ## extract limits etc. from input mask
        spacing <- attr(mask,'area')^0.5 * 100
        bb <- attr(mask, 'boundingbox')
        maxx <- max(bb$x)
        minx <- min(bb$x)
        maxy <- max(bb$y)
        miny <- min(bb$y)
        nx <- round( (maxx - minx) / spacing )
        ny <- round( (maxy - miny) / spacing )
        n <- nx * ny

        ## create rasterLayer
##        if (!require(raster))
##            stop ("unable to load raster package")
        layer <- raster(nrows = ny, ncols = nx, xmn = minx, xmx = maxx,
                                ymn = miny, ymx = maxy)

        ## A. Percolation map generation
        values(layer) <- rep(0, n)
        values(layer)[sample.int(n, round(n*p))] <- 1   ## as close as possible
        if (plt) plot(layer, useRaster = FALSE)

        ## B. Cluster identification (single-linkage clustering of adjoining pixels)
        clumped <- clump(layer, directions = directions, gaps = FALSE)

        ## C. Cluster type assignment
        ncluster <- max(values(clumped), na.rm = TRUE)
        types <- factor(c(0,1))          # nonhabitat, habitat
        numTypes <- as.numeric(types)    # 1,2
        clustertype <- sample(numTypes, ncluster, replace = TRUE, prob = c(1-A,A))
        values(clumped) <- clustertype[values(clumped)]
        if (plt) plot(clumped, useRaster = FALSE)

        ## D. Filling in image
        cellsUnassigned <- (1:n)[is.na(values(clumped))]
        cellsAssigned <- (1:n)[!is.na(values(clumped))]
        tempadj <- adjacent (clumped, cells = cellsUnassigned, target = cellsAssigned,
                             directions = 8)
        tempadj <- split(tempadj[,2], tempadj[,1])
        fillinType <- function (adjcells) {
            type <- values(clumped)[adjcells]
            type <- factor(type)
            freq <- tabulate(type)
            result <- as.numeric(levels(type)[freq == max(freq)])
            if (length(result)>1)
                result <- sample (result, 1)
            result
        }
        # cells with typed neighbours
        filled <- sapply(tempadj, fillinType)
        filledCells <- as.numeric(names(filled))
        values(clumped)[filledCells] <- filled
        # cells with no typed neighbours
        notfilledCells <- cellsUnassigned[!(cellsUnassigned %in% filledCells)]
        randomType <- sample(numTypes, length(notfilledCells), replace = TRUE, prob = c(1-A,A))
        values(clumped)[notfilledCells] <- randomType
        if (plt) plot(clumped, useRaster = FALSE)
        values(clumped) <- values(clumped) - 1

        ## optionally filter small patches
        if (minpatch > 1) {
            reclumped <- clump(clumped, directions = directions, gaps = FALSE)
            nperclump <- table(values(reclumped))
            values(clumped)[nperclump[as.character(values(reclumped))] < minpatch] <- 0
            temp <- clumped
            values(temp) <- 1-values(temp)   ## swap 0,1
            reclumped <- clump(temp   , directions = directions, gaps = FALSE)
            nperclump <- table(values(reclumped))
            values(clumped)[nperclump[as.character(values(reclumped))] < minpatch] <- 1
        }

        # pad if necessary (sets attribute OK with padding cells FALSE)
        rectmask <- rectangularMask(mask)
        # restrict to 'habitat'; discard padding
        tempmask <- subset(rectmask, attr(rectmask, 'OK') & (values(clumped)>0))
        # optionally return mask with only 'habitat' points
        if (drop) {
            mask <- tempmask
            attr(mask, 'type') <- paste('MRC p=',p, ' A=',A, sep='')
        }
        else {
            hab <-  as.numeric(pointsInPolygon (mask, tempmask))
            if (is.null(covariates(mask))) {
                covariates(mask) <- data.frame(hab)
                names(covariates(mask)) <- covname
            }
            else
                covariates(mask)[,covname] <- hab
        }
        mask
    }
}
