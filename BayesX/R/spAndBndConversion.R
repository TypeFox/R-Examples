#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* campus *.* lmu *.* de]
## Project: BayesX
## Time-stamp: <[spAndBndConversion.R] by DSB Die 09/06/2009 21:19 (CEST)>
##
## Description:
## Convert a SpatialPolygons object from package sp to the bnd format
## required by BayesX, and vice versa.
##
## History:
## 19/02/2009   file creation
## 22/02/2009   explicit scoping of sp package functions
#####################################################################################


bnd2sp <- function(bndObject)
{
    ## check if S3 class of bndObject is "bnd"
    stopifnot(inherits(x=bndObject,
                       what="bnd"))

    ## extracts
    bndNames <- names(bndObject)
    regions <- unique(bndNames)
    bndAttributes <- attributes(bndObject)
    
    ## close all polygons (last coordinates must match first coordinates)
    bndObject <- lapply(bndObject,
                        FUN=function(polygon){
                            if(! isTRUE(identical(polygon[1, ],
                                                  polygon[nrow(polygon), ])))
                            {
                                rbind(polygon,
                                      polygon[1, ])
                            } else {
                                polygon
                            }
                        })
    
    ## set up return list
    ret <- list()
    
    ## process all unique regions
    for(id in regions)
    {
        ## which polygons belong to this region?
        idMatches <- which(id == bndNames)

        ## convert these polygons to Polygon class objects
        idPolygons <- lapply(bndObject[idMatches],
                             FUN=sp::Polygon,
                             hole=FALSE)

        ## add the Polygons object with these Polygon parts to return list
        ret[[id]] <- sp::Polygons(srl=idPolygons,
                                  ID=id)        
    }

    ## add holes of inner polygons to outer regions
    surrounding <- bndAttributes$surrounding
    whichAreInner <- which(sapply(surrounding, length) > 0L)
    for(innerInd in whichAreInner)
    {
        ## get the hole
        hole <- sp::Polygon(coords=bndObject[[innerInd]],
                            hole=TRUE)

        ## get outer polys list
        outerId <- surrounding[[innerInd]]
        outerPolys <- ret[[outerId]]@Polygons

        ## add the hole to outer polys
        outerPolys <- c(outerPolys,
                        hole)

        ## write back extended outer polys list as new Polygons object with same ID as before 
        ret[[outerId]] <- sp::Polygons(srl=outerPolys,
                                       ID=outerId)
    }
    
    ## convert list of Polygons to a SpatialPolygons object and return that
    ret <- sp::SpatialPolygons(Srl=ret)
    return(ret)
}



sp2bnd <-
    function(spObject,            # object of class SpatialPolygons (or specializations)
             regionNames=sapply(spObject@polygons, slot, "ID"), # character vector of region names
                                        # (parallel to the Polygons list in spObject)
             height2width=round(diff(sp::bbox(spObject)[2, ]) / diff(sp::bbox(spObject)[1, ]), 2),
                                        # ratio of height to width
             epsilon=sqrt(.Machine$double.eps)) # how much can two polygons differ (in maximumn
                                        # squared Euclidean distance) and still match each other?
    
{
    ## check if S4 class of spObject is "SpatialPolygons"
    stopifnot(is(object=spObject,
                 class2="SpatialPolygons"))

    ## extracts
    spObject <- sp::polygons(spObject)      # now surely a SpatialPolygons object
    spList <- spObject@polygons         # discard other slots
    nRegions <- length(spList)

    ## check if number of regions matches with the length of regionNames etc
    stopifnot(is.character(regionNames),
              identical(length(regionNames), nRegions),
              height2width > 0)

    ## set up return and holes list
    ret <- list()
    holes <- list()    

    ## number of polygons and holes already processed
    numPolysProcessed <- 0
    numHolesProcessed <- 0
    
    ## process each region
    for(regionIterator in seq_along(spList))
    {
        thisRegion <- spList[[regionIterator]]@Polygons
        
        ## process each Polygon in this region
        for(polygonObject in thisRegion)
        {
            ## if it is a hole, put it in holes, else in ret.
            ## the name is set to the region name so we know from which region this
            ## polygon stems.
            if(polygonObject@hole)
            {
                ## increment hole counter
                numHolesProcessed <- numHolesProcessed + 1

                ## and correct invariant
                holes[[numHolesProcessed]] <- sp::coordinates(polygonObject)
                names(holes)[numHolesProcessed] <- regionNames[regionIterator]
            } else {
                ## increment Polygon counter
                numPolysProcessed <- numPolysProcessed + 1

                ## and correct invariant
                ret[[numPolysProcessed]] <- sp::coordinates(polygonObject)
                names(ret)[numPolysProcessed] <- regionNames[regionIterator]
            }
        }
    }
    ## sanity check
    stopifnot(all.equal(length(ret), numPolysProcessed),
              all.equal(length(holes), numHolesProcessed))

    ## now process all holes:

    ## set up surrounding list
    surrounding <- replicate(n=numPolysProcessed,
                             character())
    
    ## use number of coordinates as hash for quicker search for the embedded region
    polyDims <- sapply(ret, nrow)
    holeDims <- sapply(holes, nrow)

    for(i in seq_along(holes))
    {
        ## hash lookup
        potentialMatchesInds <- which(holeDims[i] == polyDims)

        ## now more precise search in these potential matches
        matchFound <- FALSE
        for(j in potentialMatchesInds)
        {
            ## decide
            thisHole <- holes[[i]]
            thisRegion <- ret[[j]]

            squaredEuclideanDistances <-
                rowSums((thisHole[rev(seq_len(nrow(thisHole))), ] - thisRegion)^2) 
            doesMatch <- max(squaredEuclideanDistances) < epsilon

            ## if it matches, update the surrounding data
            ## and break out of the for loop
            if(doesMatch)
            {
                matchFound <- TRUE

                surrounding[[j]] <- names(holes)[i]
                
                ## we can proceed with the next hole:
                break
            }
        }

        ## echo a warning if a hole has no match
        if(! matchFound)
        {
            warning(simpleWarning(paste("No match found for hole in region",
                                        names(holes)[i])))
        }
    }

    ## finally collect information and return the bnd object
    ret <- structure(ret,
                     surrounding=surrounding,
                     height2width=height2width,
                     class="bnd")
    return(ret)
}

