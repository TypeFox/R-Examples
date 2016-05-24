"export.asc" <- function(x, file)
{
    ## verifications
    if (!inherits(x, "asc")) stop("Non convenient data")
    if (substr(file, nchar(file)-3, nchar(file))!=".asc")
        file<-paste(file, ".asc", sep="")

    ## Creates the file header
    file.create(file)
    zz<-file(file, "w")
    nc<-paste("ncols", "         ", nrow(x), sep="")
    nl<-paste("nrows", "         ", ncol(x), sep="")
    xll<-paste("xllcorner", "     ",
               attr(x, "xll")-attr(x, "cellsize")/2, sep="")
    yll<-paste("yllcorner", "     ",
               attr(x, "yll")-attr(x, "cellsize")/2, sep="")
    cs<-paste("cellsize", "      ", attr(x, "cellsize"), sep="")
    nas<-paste("NODATA_value", -9999, sep="  ")

    ## write to the file
    writeLines(nc, zz)
    writeLines(nl, zz)
    writeLines(xll, zz)
    writeLines(yll, zz)
    writeLines(cs, zz)
    writeLines(nas, zz)

    close(zz) ## close the connection


    ## replace the missing values, adds newlines at the end
    ## of the rows OF THE MAP (so column of the matrix)
    x[is.na(x)]<--9999
    x<-x[,ncol(x):1]
    x<-rbind(x, rep("\n", ncol(x)))

    ## ... and sinks to the file
    sink(file, append=TRUE)
    cat(x)
    sink()

  }

