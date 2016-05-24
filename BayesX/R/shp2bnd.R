shp2bnd <- function(shpname, regionnames, check.is.in = TRUE)
{   
    ## safe coercions ...
    shpname <- as.character(shpname)
    regionnames <- as.character(regionnames)
    check.is.in <- as.logical(check.is.in)

    ## ... and checks
    stopifnot(identical(length(shpname), 1L),
              length(regionnames) >= 1L,
              identical(length(check.is.in), 1L))
    
    ## now read the shapefile information
    shp <- shapefiles::read.shapefile(shpname)
    dbf <- shapefiles::read.dbf(paste(shpname,".dbf",sep=""))

    ## extract names of the regions:
    regionnames <-
        ## if there is only one string, then we assume it is the variable name
        ## in the database of the shape file
        if(identical(length(regionnames), 1L))
        {
            ## so get that variable
            as.character(dbf$dbf[regionnames][[1L]])
        } else {
            ## we stay at the given names
            regionnames
        }
    
    ## delete commas in names of the regions
    newRegionnames <- gsub(pattern=",",
                           replacement="",
                           x=regionnames,
                           fixed=TRUE)

    ## check if commas have been deleted and if so, issue a warning
    if(! identical(regionnames, newRegionnames))
        warning(simpleWarning("commas in names of the regions have been deleted"))

    ## overwrite old region names with new region names
    regionnames <- newRegionnames
    rm(newRegionnames)

    ## split data into closed polygons.
    ## we need:
    
    ## storage for the new names and the new polygons
    polyList <- list()
    ## and the corresponding original indexes
    originalRegionNumber <- integer()

    cat("Reading map ...")
    
    ## now process all original regions
    for(i in seq_along(regionnames))
    {
        ## this is temporary storage (originally a data.frame with X and Y columns):
        temppoly <- as.matrix(shp$shp$shp[[i]]$points)
        dimnames(temppoly) <- NULL

        ## as long as there are still points to be processed
        while((nPoints <- nrow(temppoly)) > 0L)
        {
            ## where does the first point occur in the data at the second time?
            endIndex <- which((temppoly[- 1L, 1] == temppoly[1L, 1]) &
                              (temppoly[- 1L, 2] == temppoly[1L, 2])) + 1L

            ## take the first next occurrence, or the last point if the polygon is not closed
            endIndex <- 
                if(length(endIndex) > 0L)
                {
                    endIndex[1L]
                } else {
                    nPoints
                }           

            ## the range of this polygon
            polyRange <- 1L:endIndex

            ## this was index i
            originalRegionNumber <- c(originalRegionNumber,
                                      i)
            
            ## save the polygon
            polyList <- c(polyList,
                          list(temppoly[polyRange, ]))
            ## list is necessary so that c(list(), data.frame(..)) is a one-element list,
            ## and not a list with the variables of the data.frame as elements

            ## and delete this part from temporary storage
            temppoly <- temppoly[- polyRange, ]
        }
    }

    cat(" finished\n")
    
    ## so how many polygons do we have now?
    nPolys <- length(polyList)

    cat("Note: map consists originally of", nPolys, "polygons\n")

    ## here is the parallel list of the surrounding region names of single polygons
    surrounding <- replicate(n=nPolys,
                             expr=character()) ## until now no region names anywhere!
    
    ## check for polygons contained in another polygon?
    if(check.is.in)
    {      
        ## get dimensions of all polygons
        dims <- sapply(polyList, nrow)

        ## save here which polygons should be removed, because they are boundaries
        ## to polygons lying inside
        rmcheck <- logical(nPolys)

        ## save here the indexes of the polygons which have already been matched/processed.
        ## these must not be processed again!
        whichWereProcessed <- integer()

        ## process each polygon i
        for(i in seq_len(nPolys))
        {
            ## if we had processed this already
            if(i %in% whichWereProcessed)
            {
                ## go on to the next polygon
                next
            } else {
                ## add i to processed ones
                whichWereProcessed <- union(whichWereProcessed,
                                            i)
            }
            
            ## which polygons have same number of points as the current?
            sameDimsInds <- setdiff(which(dims == dims[i]),
                                    whichWereProcessed) ## but without the already processed ones     

            ## process all polygons j with same dims as polygon i
            ## (this works as a hash)
            for(j in sameDimsInds)
            {
                ## compute squared distance of polygon_i and reversed polygon_j
                reverseInds <- dims[i]:1L
                squaredDistance <- sum( (polyList[[i]] -
                                         polyList[[j]][reverseInds, ])^2 )

                ## if it is small enough
                if(squaredDistance < 1e-5)
                {                       
                    ## find out which is the outer one
                    outer <- inner <- 0L
                    if(.ringDirxy(polyList[[j]]) < 0)
                    {
                        outer <- j
                        inner <- i
                    } else {
                        outer <- i
                        inner <- j
                    }

                    ## remove the outer polygon
                    rmcheck[outer] <- TRUE

                    ## and add the information in which region it is lying
                    ## (each polygon can only lie in 1 other region, of course...)
                    surrounding[[inner]] <- regionnames[originalRegionNumber[outer]]
                }

                ## we have processed j
                whichWereProcessed <- union(whichWereProcessed,
                                            j)
            }            
        }

        ## we have processed all polygons, and can remove the unnecessary ones
        polyList <- polyList[! rmcheck]
        originalRegionNumber <- originalRegionNumber[! rmcheck]
        surrounding <- surrounding[! rmcheck]

        cat("Note: After removing unnecessary surrounding polygons, the map consists of",
            length(polyList), "polygons\n")          
    }    

    ## add the original region names to the polygons list as names
    names(polyList) <- regionnames[originalRegionNumber]

    ## the new unique regions
    regions <- unique(names(polyList))
    cat("Note: map consists of", length(regions), "regions\n")
    
    ## compute relation of height to width (for plotting etc)
    minima <- sapply(polyList, function(x){apply(x,2,min)})
    maxima <- sapply(polyList, function(x){apply(x,2,max)})
    
    minimum <- apply(minima,1,min)
    maximum <- apply(maxima,1,max)
    
    x.range <- maximum[1] - minimum[1]
    y.range <- maximum[2] - minimum[2]
    
    height2width <- round(y.range / x.range, digits=2)
    
    ## now return the bnd object
    return(structure(polyList,
                     class="bnd",
                     height2width=height2width,
                     surrounding=surrounding,
                     regions=regions))
}
