bnd2gra <- function (map)
{
    ## check class of argument
    if (! inherits(map, "bnd"))
        stop("argument 'map' is not an object of class 'bnd'")

    ## extract (unique) regions from the polygons
    regions <- unique(names(map))
    nRegions <- length(regions)

    ## initialize return matrix
    retMatrix <- matrix(data=0,
                        nrow=nRegions, ncol=nRegions,
                        dimnames=list(regions, regions))

    ## helper function:
    pointsMatrixToPointsCharVector <- function(pointsMatrix)
    {
        paste(pointsMatrix[, 1L],
              pointsMatrix[, 2L],
              sep="&")
    }

    ## check for neighboring polygons:
    ## process only upper triangular part of the matrix, because it is symmetric.

    ## counter for already processed pairs
    nPairsProcessed <- 0L

    ## how many pairs need to be processed?
    nPairs <- (nRegions * (nRegions - 1L)) / 2L

    ## and at which iterations should a progress message be printed?
    nProgressIterations <- 10L
    progressIterations <- seq(from=0L,
                              to=nPairs,
                              length=nProgressIterations + 1L)[- 1L]

    cat("Start neighbor search ...\n")
    
    pointsMatrix <- NULL

    ## i is the row region number
    for (i in seq_len(nRegions - 1L))
    {
        ## which polygons belong to region number i?
        polyIndsRegion.i <- which(names(map) == regions[i])

        ## make a vector of all points (x, y) in the format "x&y" which belong to this region
        if(length(polyIndsRegion.i)>1)
          {
          pointsMatrix <- map[polyIndsRegion.i[1]]
          for(k in 2:length(polyIndsRegion.i))
            pointsMatrix[[1]] <- rbind(pointsMatrix[[1]], map[polyIndsRegion.i[k]][[1]])
          }
        else 
          {
          pointsMatrix <- map[polyIndsRegion.i]
          }
        pointsRegion.i <- sapply(X=pointsMatrix,
                                 FUN=pointsMatrixToPointsCharVector)

        ## j is the column region number
        for (j in (i + 1):nRegions)
        {
            ## which polygons belong to region number j?
            polyIndsRegion.j <- which(names(map) == regions[j])
            
            ## make a vector of all points (x, y) in the format "x&y" which belong to this region
            if(length(polyIndsRegion.j)>1)
              {
              pointsMatrix <- map[polyIndsRegion.j[1]]
              for(k in 2:length(polyIndsRegion.j))
                pointsMatrix[[1]] <- rbind(pointsMatrix[[1]], map[polyIndsRegion.j[k]][[1]])
              }
            else 
              {
              pointsMatrix <- map[polyIndsRegion.j]
              }
            pointsRegion.j <- sapply(X=pointsMatrix,
                                     FUN=pointsMatrixToPointsCharVector)

            ## now decide if region i and j share at least 2 common points in their polygons
            if (sum(pointsRegion.i %in% pointsRegion.j) >= 2L)
            {
                ## then they are neighbors!
                retMatrix[i, j] <-
                    retMatrix[j, i] <- - 1L
            }

            ## increment counter
            nPairsProcessed <- nPairsProcessed + 1L
            
            ## echo progress?
            if (nPairsProcessed %in% progressIterations)
            {
                ## echo percentage of processed pairs
                percentageProcessed <- floor((nPairsProcessed * 100L) / nPairs)
                cat(paste("progress: ", percentageProcessed, "%\n",
                          sep = ""))
            }            
        }
    }

    cat("Neighbor search finished.\n")

    ## add is.in relations
    surrounding <- attr(map, "surrounding")
    whichPolygonsAreInner <- which(sapply(surrounding, length) > 0L)   
    for(innerInd in whichPolygonsAreInner)
    {
        innerRegion <- names(map)[innerInd]
        outerRegion <- surrounding[[innerInd]]
        
        retMatrix[innerRegion, outerRegion] <-
            retMatrix[outerRegion, innerRegion] <- - 1L
    }

    ## on the diagonal, there are the negative row sums
    diag(retMatrix) <- - rowSums(retMatrix)

    ## finally return the matrix as a graph object
    class(retMatrix) <- "gra"
    return(retMatrix)
}
