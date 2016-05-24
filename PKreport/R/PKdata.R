#############################################################################################
## File: PKdata.R
## Author: Xiaoyong Sun
## Date: 10/20/2009
## Goal: read in PK data; file config
## Notes:
##      -
#############################################################################################


PKdata <- function(data, match.term=NULL)
{
    if (missing(data)) stop("Data is required!")
    if (length(data)==0 || nrow(data)==0) stop("Data is Not available!")
    
    ## read in data
    if (length(colnames(data))==0) stop("Data column does NOT have names!")

    # check missing value
    sapply(1:ncol(data), function(i)
          {
              if (all(is.na(data[,i]))) stop(paste("\nData column ", i, " are all NA values!", sep=""))
              if (any(is.na(data[,i]))) warning(paste("\nData column ", i, " has missing values!", sep=""))
          })
    
    ## match term
    if (is.null(match.term)) stop("Please input config list!")
    mt <- unlist(match.term)
    PK.match <- match(mt, colnames(data))
    if(length(PK.match[is.na(PK.match)]) > 0) stop(paste(dQuote(mt[is.na(PK.match)]) , "in config list do NOT match data!\n", sep=" "))

    ## make sure WRES, RES, PRED, IPRE, DV, TIME are only one item ??
    # check ID, DV, TIME
    if ( (length(match.term$ID)!=1) || (length(match.term$DV)!=1)  || (length(match.term$TIME)!=1) )
    stop("Please make sure ID, DV and TIME are input with only ONE variable!")

    if ( (length(match.term$RES)!=1) || (length(match.term$WRES)!=1) )
    stop("Please make sure RES and WRES are input with only ONE variable!")

    if ( (length(match.term$PRED)!=1) || (length(match.term$IPRE)!=1) )
    stop("Please make sure PRED and IPRE are input with only ONE variable!")
    
    .pkplot$setTerm(match.term)
    .pkplot$setPKData(data)
                                  
    cat("Data is read successfully.\n")

}

 # NOTE: "general.list" should be "global.list"
PKconfig <- function(general.list, hist.list, scatter.list)
{
    if (any(!(general.list$save.format %in% c("jpeg", "bmp", "png", "tiff", "win.metafile"))))
    {
        stop("The save format is NOT supported! jpeg, bmp, png, tiff, win.metafile are supported!")
    }

    # check general.list
    if (is.null(general.list$save.format) || !("png" %in% general.list$save.format))
    {
        general.list$save.format <- c("png", unique(general.list$save.format))
    }
    else
    {
        png.ind <- which("png" == general.list$save.format)
        general.list$save.format <- c("png", general.list$save.format[-png.ind])
    }


    ## general term setup
    sapply(names(general.list), function(i) .pkplot$setGlobalConfig(i,general.list[[i]]))
    
    ## graph global term
    ## - for lattice
    lattice.global <- c("col", "span", "type", "layout")
    ggplot.global <- c("col", "span")

    ## hist - lattice setup config term
    lattice.list <- hist.list[names(hist.list) %in% lattice.global]
    .pkplot$setHistGraph(lattice.list, "lattice")

    ## hist - ggplot setup config term
    ggplot.list <- hist.list[names(hist.list) %in% ggplot.global]
    ggplot.list$geom <- c("histogram")
    .pkplot$setHistGraph(ggplot.list, "ggplot")
    
    ## scatter - lattice setup config term
    lattice.list <- scatter.list[names(scatter.list) %in% lattice.global]
    .pkplot$setScatterGraph(lattice.list, "lattice")

    ## scatter - ggplot setup config term
    ggplot.list <- scatter.list[names(scatter.list) %in% ggplot.global]
    ## ggplot type setup || layout setup too
    if (!is.null(scatter.list$type))
    {
        ggplot.list$geom <- switch(paste(scatter.list$type, collapse=""),
                                    p = c("point"),
                                    l = c("line"),
                                    psmooth = c("point", "smooth"),
                                    lsmooth = c("line", "smooth"))
    }
    else
    {
        warning("ggplot package does not have matching type at this time!")
    }

    .pkplot$setScatterGraph(ggplot.list, "ggplot")

    others.list <- hist.list[!(names(hist.list) %in% lattice.global)]
    .pkplot$setHistGraph(others.list, "others")
    others.list <- scatter.list[!(names(scatter.list) %in% lattice.global)]
    .pkplot$setScatterGraph(others.list, "others")

    invisible(NULL)
    
}