setMethod(
    f = "stratify",
    signature = signature(
        object = "SpatialPixels"
    ),
    definition = function(object, nStrata, priorPoints = NULL, maxIterations = 1000L, nTry = 1L,
        equalArea = FALSE, verbose = getOption("verbose")) {

        # coerce object to class "SpatialPixels" (needed for descendants of "SpatialPixels" like "SpatialGridDataFrame")
        #gridded(object) <- FALSE
        #object <- suppressWarnings(as(as(object, "SpatialPoints"), "SpatialPixels"))

        # check 'nStrata' argument
        if (length(nStrata) != 1L) {
            stop("'nStrata' should be an integer vector of length 1 (i.e., a scalar)", call. = FALSE)
        }
        if (nStrata < 1L) {
            stop("'nStrata' should be a strictly positive integer", call. = FALSE)
        }

        # check if prior points have been specified
        hasPriorPoints <- !is(priorPoints, "NULL")

        # check if prior points is an instance of class "SpatialPoints"
        if (hasPriorPoints) {
            if (class(priorPoints) != "SpatialPoints") {
                stop("'priorPoints' has to be an instance of class \"SpatialPoints\"", call. = FALSE)
            }
            if (nrow(coordinates(priorPoints)) == 0L) { # rare, but not impossible
                stop("'priorPoints' does not contain coordinates", call. = FALSE)
            }
        }

        # check if a projection have been used
        externalProjection <- proj4string(object)
        hasProjection <- !is.na(externalProjection)

        # convert current projection to lat-long WGS84
        if (hasProjection) {
            if (equalArea) {
                stop("strata of equal area in combination with map\nprojections are currently not supported ", call. = FALSE)
            }
            if (suppressWarnings(require(rgdal))) {
                internalProjection <- CRS("+proj=longlat +ellps=WGS84")
                object <- spTransform(object, internalProjection)
                suppressWarnings(gridded(object) <- TRUE)
            } else {
                stop("You need package 'rgdal' to handle projection attributes\nPlease install this package first", call. = FALSE )
            }
        }

        # extract coordinates of cell centers
        sObject <- coordinates(object)

        # get and check the number of grid cells
        nCells <- nrow(sObject)
        if (nCells < nStrata) {
            stop("'object' has insufficient grid cells", call. = FALSE)
        }

        # Initialize the probability of selecting a grid cell as initial cluster center.
        # These probabilities will all be equal, except for cells containing prior points.
        # The probability of these cells will be set to zero.
        # (only needed for the transfer algorithm, not for the swopping algorithm)
        if (!equalArea) {
            cellSelectionProbability <- rep(x = 1, times = nCells)
        }

        # check and process prior points
        if (hasPriorPoints) {

            # check 'equalArea' argument
            if (equalArea) {
                stop("'equalAreas' should be FALSE in case of prior points", call. = FALSE)
            }

            # check and transform projection
            if (hasProjection) {
                if (!identical(externalProjection, proj4string(priorPoints))) {
                    stop("projections of 'object' and 'priorPoints' don't match", call. = FALSE)
                }
                priorPoints <- spTransform(priorPoints, CRS("+proj=longlat +ellps=WGS84"))
            }

            # Remove all prior points outside the target universe.
            # These points may lead to empty clusters.
            # Also remove prior points that are nearly coinciding.
            # Nearly coinciding points will result in empty clusters.
            # Points are said to be nearly coinciding if they are in the same grid cell.
            cellIdPriorPoint <- priorPoints %over% geometry(object)
            isInTargetUniverse <- !is.na(cellIdPriorPoint)
            isCoinciding <- duplicated(cellIdPriorPoint)
            keep <- isInTargetUniverse & !isCoinciding
            priorPoints <- priorPoints[keep, , drop = FALSE]
            cellIdPriorPoint <- cellIdPriorPoint[keep]
            if (any(!keep)) {
                if (any(!isInTargetUniverse)) {
                    warning(sum(!isInTargetUniverse), " location(s) outside the target universe ",
                        "(as defined by 'object') have been found\n",
                        "These locations have been removed. ", sum(isInTargetUniverse),
                        " location(s) have been retained", call. = FALSE)
                }
                if (any(isCoinciding)) {
                    warning("(Nearly) coinciding points have been removed. ", call. = FALSE)
                }
            }

            # The initial set of centroids consists of prior points and a set of randomly selected points
            # Make sure that the latter don't coincide with the latter by setting the selection probability
            # of grid cells containing prior points to zero.
            cellSelectionProbability[cellIdPriorPoint] <- 0

            # extract coordinates
            sPriorPoints <- coordinates(priorPoints)

            # make sure that nStrata is always greater than or equal
            # to the number of prior points
            nPriorPoints <- nrow(sPriorPoints)
            if (nStrata < nPriorPoints) {
                stop("'nStrata' should be greater than or equal to ",
                    "the number of 'priorPoints'", call. = FALSE)
            }
        } else {
            nPriorPoints <- 0L
            sPriorPoints <- NULL
        }

        # set number of free points (i.e., cluster centers to optimize)
        nFreePoints <- nStrata - nPriorPoints

        # create instances of cell centers
        cellCenters <- J(
            class = "partition/LocationFactory",
            method = ifelse(hasProjection, "getLatLongInstance", "getEastingNorthingInstance"),
            .jarray(as(sObject[, 1], "double")),
            .jarray(as(sObject[, 2], "double"))
        )

        # initialize objective function value of optimal configuration
        minimumObjectiveFunctionValue <- Inf

        # start loop with different initial configurations
        for (i in seq_len(nTry)) {

            # show progress
            if (verbose) {                
                cat(format(Sys.time()), "| optimizing configuration", i, "\n")
            }

            # create a Java reference to a subclass of class Stratification
            if (equalArea) {
                # create an instance of class "CompactSpatialPartitionSwop"
                p <- new(
                    J("partition/CompactSpatialPartitionSwop"),
                    cellCenters,
                    .jarray(as(sample(x = rep(x = 0:(nStrata - 1), length = nCells)), "integer"))
                )
            } else {

                # select initial cluster centers (random sampling of cells without replacement)
                # cells that contain a prior point will be excluded to prevent collocation
                cellId <- sample(x = seq_len(nCells), size = nFreePoints, replace = FALSE, prob = cellSelectionProbability)
                sClusterCenters <- rbind(sPriorPoints, sObject[cellId, ])

                # create instances of cluster centers
                clusterCenters <- J(
                    class = "partition/LocationFactory",
                    method = ifelse(hasProjection, "getLatLongInstance", "getEastingNorthingInstance"),
                    .jarray(as(sClusterCenters[, 1], "double")), #.jarray forces scalar arguments to arrays of length 1
                    .jarray(as(sClusterCenters[, 2], "double"))
                )

                # create an instance of class "CompactSpatialPartitionTransfer"
                p <- new(
                     J("partition/CompactSpatialPartitionTransfer"),
                     cellCenters,
                     clusterCenters,
                     .jarray(rep(x = c(TRUE, FALSE), times = c(nPriorPoints, nFreePoints)))
                )
            }

            # set maximum number of iterations
            p$setMaximumNumberOfIterations(as(maxIterations, "integer"))

            # optimize partition
            p$optimize()

            # retrieve objective function value
            objectiveFunctionValue <- p$getObjectiveFunctionValue()

            # check if current stratification is more optimal than previous stratification(s)
            if (objectiveFunctionValue < minimumObjectiveFunctionValue) {
                hasConverged <- p$hasConverged()
                minimumObjectiveFunctionValue <- objectiveFunctionValue
                clusterId <- as(p$getClusterId(), "integer")
                centroids <- p$getCentroids()
                if (hasProjection) {
                    sCentroids <- data.frame(
                        s1 = sapply(X = centroids, FUN = function(x) {x$getLongitude()}),
                        s2 = sapply(X = centroids, FUN = function(x) {x$getLatitude()})
                    )
                } else {
                    sCentroids <- data.frame(
                        s1 = sapply(X = centroids, FUN = function(x) {x$getEasting()}),
                        s2 = sapply(X = centroids, FUN = function(x) {x$getNorthing()})
                    )
                }
            }
            
            # show progress
            if (verbose) {                
                cat(format(Sys.time()), "|    current objective function value:", objectiveFunctionValue, "\n")
                cat(format(Sys.time()), "|    minimum objective function value:", minimumObjectiveFunctionValue, "\n")
            }
        }

        # check convergence
        if (!hasConverged) {
            warning("Convergence has not been reached after ",
                maxIterations, " interations and ", nTry,
                " attempts", call. = FALSE)
        }
        
        # promote sCentroids to class "SpatialPoints"
        colnames(sCentroids) <- colnames(sObject)
        coordinates(sCentroids) <- colnames(sCentroids)
        centroids <- sCentroids

        # restore original projection
        if (hasProjection) {
            object <- spTransform(object, CRS(externalProjection))
            suppressWarnings(gridded(object) <- TRUE)
            proj4string(centroids) <- internalProjection
            centroids <- spTransform(centroids, CRS(externalProjection))
            if (hasPriorPoints) {
                priorPoints <-  spTransform(priorPoints, CRS(externalProjection))
            }
        }

        # create an instance of a subclass of  "Stratification"
        if (equalArea) {
            stratification <- new(
                Class = "CompactStratificationEqualArea",
                cells = object,
                stratumId = clusterId,
                centroids = centroids,
                mssd = minimumObjectiveFunctionValue
            )
        } else {
            if (hasPriorPoints) {
                stratification <- new(
                    Class = "CompactStratificationPriorPoints",
                    cells = object,
                    stratumId = clusterId,
                    centroids = centroids,
                    priorPoints = priorPoints,
                    mssd = minimumObjectiveFunctionValue
                )
            } else {
                stratification <- new(
                    Class = "CompactStratification",
                    cells = object,
                    stratumId = clusterId,
                    centroids = centroids,
                    mssd = minimumObjectiveFunctionValue
                )
            }
        }

        # return stratification
        stratification
    }
)
