setMethod(
    f = "spsample",
    signature = signature(
        x = "CompactStratificationEqualArea",
        n = "numeric",
        type = "character"
    ),
    definition = function(x, n, type, ...) {
        samplingPattern <- spsample(x = x, n = n, ...)
        type <- match.arg(arg = type, choices = "composite", several.ok = FALSE)
        if (type == "composite") {
            numberOfStrata <- getNumberOfStrata(x)
            samplingPattern <- new(
                Class = "SamplingPatternRandomComposite",
                sample = as(samplingPattern, "SpatialPoints"),
                composite = rep(x = seq_len(n), times = numberOfStrata)
            )
        }
        samplingPattern
    }
)



setMethod(
    f = "spsample",
    signature = signature(
        x = "CompactStratification",
        n = "missing",
        type = "missing"
    ),
    definition = function(x, ...) {
        
        # extract the number of strata
        nStrata <- getNumberOfStrata(x)
        
        # extract centroids
        centroids <- getCentroid(x)
        sCentroids <- coordinates(centroids)
        sNames <- colnames(sCentroids)
        
        # Assign centroids outside the target universe to the nearest
        # cell within the target universe. To simplify things, the Euclidean
        # distance will used until a better solution has been found for
        # handling these kinds of centroids
        isOutsideTargetUniverse <- is.na(centroids %over% geometry(x@cells))
        if (any(isOutsideTargetUniverse)) {
            sCells <- t(coordinates(x@cells))
            for (i in which(isOutsideTargetUniverse)) {
                squaredDistance <- colSums((sCells - sCentroids[i, ])^2)
                j <- which.min(squaredDistance)
                sCentroids[i, ] <- sCells[, j]
            }
            centroids <- as.data.frame(sCentroids)
            coordinates(centroids) <- sNames
        }
        
        # return an object of class "SamplingPattern"
        new(
            Class = "SamplingPatternCentroids",
            sample = centroids
        )
    }
)



setMethod(
    f = "spsample",
    signature = signature(
        x = "CompactStratification",
        n = "numeric",
        type = "missing"
    ),
    definition = function(x, n, ...) {
        
        # extract cell size
        cellSize <- getCellSize(x)
        
        # extract coordinates
        s <- coordinates(x@cells)
        
        # randomly select 'n' cells per stratum _with_ replacement
        cellId <- tapply(
            X = seq_len(nrow(s)),
            INDEX = x@stratumId,
            FUN = function(cell) {
                if (length(cell) > 1) {
                    return(sample(x = cell, size = n, replace = TRUE))
                } else { # in case only one cell is available
                    return(rep(x = cell, times = n))
                }
            }
        )
        cellId <- unlist(x = cellId, use.names = FALSE)
        s <- s[cellId, , drop = FALSE]
        
        # randomly select one location in each selected cell
        U <- runif(n = 2 * nrow(s), min = -0.5, max = 0.5)
        s0 <- matrix(
            data =  U * cellSize, # so cells may be rectangular
            nrow = nrow(s),
            ncol = ncol(s),
            byrow = TRUE
        )
        s <- s + s0
        
        # return result as an instance of class "SamplingPatternRandomSamplingUnits"
        new(
            Class = "SamplingPatternRandomSamplingUnits",
            sample = SpatialPoints(coords = s)
        )
    }
)



setMethod(
    f = "spsample",
    signature = signature(
        x = "CompactStratificationPriorPoints",
        n = "missing",
        type = "missing"
    ),
    definition = function(x, ...) {
        
        # get centroids
        centroids <- getCentroid(x)
        priorPoints <- x@priorPoints
        
        # get coordinates of centroids
        sCentroids <- coordinates(centroids)
        sPriorPoints <- coordinates(priorPoints)
        
        # get number of centroids
        nCentroids <- nrow(sCentroids)
        nPriorPoints <- nrow(sPriorPoints)
        nFreeCentroids <- nCentroids - nPriorPoints
        
        # update centroids
        sCentroids[seq_len(nPriorPoints), ] <- sPriorPoints
        
        # return object of class "SamplingPatternPriorPoints"
        new(
            Class = "SamplingPatternPriorPoints",
            sample = SpatialPoints(coords = sCentroids),
            isPriorPoint = as(rep(x = c(TRUE, FALSE), times = c(nPriorPoints, nFreeCentroids)), "logical")
        )
    }
)
