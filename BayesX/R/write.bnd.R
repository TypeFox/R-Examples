write.bnd <- function(map, file, replace=FALSE)
{
    if(! inherits(map,"bnd"))
        stop("argument 'map' is not an object of class 'bnd'")

    ## coercions, checks
    replace <- as.logical(replace)
    file <- as.character(file)
    stopifnot(identical(length(file), 1L),
              identical(length(replace), 1L))    

    ## check whether the file exists
    if(file.exists(file))
    {
        if(replace)
        {
            removeSucceeded <- file.remove(file)
            if(! removeSucceeded)
            {
                stop("file exists, but could not be removed")
            }
        } else {
            stop("specified file already exists")
        }        
    }

    myQuote <- function(string)
    {
        return(paste("\"", string, "\"",
                     sep=""))
    }
    
    ## names of the belonging regions
    belongingRegions <- names(map)
    
    ## no. of polygons
    nPolygons <- length(map)

    ## the surrounding list
    surrounding <- attr(map, "surrounding")
    
    for(i in seq_len(nPolygons))
    {
        dat <- map[[i]]
        dat <- dat[complete.cases(dat), ]
        
        temp <- paste(myQuote(belongingRegions[i]),
                      nrow(dat),
                      sep=",")
        write(temp, file, append=TRUE)
        
        if(length(outerRegionName <- surrounding[[i]]))
        {
            con <- paste("is.in",
                         myQuote(outerRegionName),
                         sep=",")
            write(con, file, append=TRUE)
        }
        write.table(dat, file, append=TRUE,
                    col.names=FALSE, row.names=FALSE,
                    sep=",", quote=FALSE)
    }

    return(invisible())
}

