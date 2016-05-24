read.bnd <- function(file, sorted=FALSE)
{
    ## list with columns from .bnd file
    data.raw <- scan(file,
                     what = list("", ""),
                     sep = ",",
                     quote = "")

    ## disable warnings and surely revert to old settings afterwards
    oldOptions <- options(warn = -1)
    on.exit(options(oldOptions))    

    ## source of warnings: NAs in region names after is.in
    data.numeric <- lapply(data.raw,
                           as.numeric)

    ## revert now to old settings to show other warnings
    options(oldOptions)

    ## helper function
    unquote <- function(string)
    {
        return(gsub(pattern="\"",
                    replacement="",
                    x=string))
    }
    
    ## where do we have is.in's?
    whereIsIn <- which(data.raw[[1]] == "is.in") 

    ## so the surrounding names are
    surroundingNames <- unquote(data.raw[[2]][whereIsIn])

    ## so where do we have the region names of the polygons?
    whereRegionNames <- setdiff(which(is.na(data.numeric[[1]])),
                                whereIsIn)

    ## and we know how many polygons there are now
    nPolygons <- length(whereRegionNames)
    cat("Note: map consists of", nPolygons, "polygons\n")

    ## extract the region names
    belongingRegions <- unquote(data.raw[[1]][whereRegionNames])

    ## so we can already setup the regions vector
    regions <- unique(belongingRegions)
    cat("Note: map consists of", length(regions), "regions\n")
    
    ## extract the lengths of the polygons
    polyLengths <- data.numeric[[2]][whereRegionNames]

    ## find out which polygon number the is.in information belongs to
    enclosedPolygonsInds <- findInterval(x=whereIsIn,
                                         vec=whereRegionNames)

    ## so we can already setup the surrounding list
    surrounding <- replicate(n=nPolygons,
                             expr=character())
    for(i in seq_along(enclosedPolygonsInds))
    {
        surrounding[[enclosedPolygonsInds[i]]] <- surroundingNames[i]
    }
    
    ## now we can cbind and delete all not-point-stuff from the numeric data
    data.numeric <- cbind(data.numeric[[1]],
                          data.numeric[[2]])
    data.numeric <- na.omit(data.numeric)
    
    ## start processing the single polygons
    cat("Reading map ...")

    ## set up the list with the correct names
    map <- vector(mode="list", length=nPolygons)
    names(map) <- belongingRegions

    ## these are the start indices for the polygons
    ## (plus one index past the last polygon)
    startInds <- cumsum(c(1, polyLengths))
    
    ## so filling the list with the points is easy now
    for(k in seq_along(map)) {
        map[[k]] <- data.numeric[startInds[k]:(startInds[k+1] - 1), ]
        if(sum(map[[k]][1,] == map[[k]][polyLengths[k],]) != 2)
           warning(paste("Note: First and last point of polygon ",k," (region ",names(map)[k],") are not identical", sep=""), call. = FALSE)
    }

    ## processing finished
    cat(" finished\n")
    rm(data.numeric)
    
    ## sort the map?
    if(sorted){
        ## try conversion of region names to numbers
        numericNames <- as.numeric(names(map))

        ## decide which ordering to apply:
        ## are there some characters?
        newOrder <- 
            if(any(is.na(numericNames)))
            {
                cat("Note: regions sorted by name\n")
                order(names(map))
                
            } else {
                cat("Note: regions sorted by number\n")
                order(numericNames)                
            }

        ## order everything 
        map <- map[newOrder]
        surrounding <- surrounding[newOrder]
    }

    ## determine aspect ratio
    minima <- sapply(map, function(x){apply(x,2,min)})
    maxima <- sapply(map, function(x){apply(x,2,max)})
    
    minimum <- apply(minima,1,min)
    maximum <- apply(maxima,1,max)
    
    x.range <- maximum[1] - minimum[1]
    y.range <- maximum[2] - minimum[2]

    ## return the bnd object
    rval <- structure(map, class = "bnd",
      surrounding = surrounding, regions = regions)
    attr(rval, "asp") <- (y.range / x.range) / cos((mean(c(maximum[2] - minimum[2])) * pi) / 180)

    rval
}

