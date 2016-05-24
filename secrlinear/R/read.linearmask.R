############################################################################################
## package 'secrlinear'
## read.linearmask.R
## last changed 2014-08-30; 2014-10-26 graph attribute;
## 2014-10-31 optional read from shapefile
## 2014-11-03 make.linearmask (called by read.linearmask and rbind.linearmask)
############################################################################################

make.linearmask <- function (SLDF, spacing, spacingfactor, graph, cleanskips)  {
    ## for bounding box...
    tmp <- lapply(coordinates(SLDF), function(x) do.call("rbind", x))
    tmp <- do.call(rbind, tmp)
    xyl <- lapply(as.data.frame(tmp), range)
    names(xyl) <- c('x','y')

    ## discretize line
    maskSPDF <- sample.line (SLDF, spacing)   ## SPDF
    if (is.null(maskSPDF)) {
        mask <- data.frame(x=numeric(0), y=numeric(0))
        covariates(mask) <- data.frame(LineID = numeric(0))
    }
    else {
         mask <- data.frame(coordinates(maskSPDF))         ## dataframe
         names(mask) <- c('x', 'y')
    }
    attr(mask, 'SLDF') <- SLDF
    attr(mask,'boundingbox') <- do.call(expand.grid, xyl)[c(1,2,4,3),]
    attr(mask,'type')    <- 'user'
    attr(mask,'spacing') <- spacing
    attr(mask,'spacingfactor') <- spacingfactor
    class(mask) <- c('linearmask', 'mask', 'data.frame')

    if (nrow(mask) > 0) {
        ## covariates
        df <- data.frame(maskSPDF)
        covariates(mask) <- df
        
        ## construct graph
        if (graph) {
            attr(mask, 'graph') <- asgraph(mask)
        }
        
        ## remove termini etc.
        OK <- !(names(df) %in% c( "coords.x1", "coords.x2", "x", "y", "Terminal"))
        terminal <- df$Terminal
        mask <- mask[!terminal,]
        covariates(mask) <- covariates(mask)[!terminal,OK]
        attr(mask,'meanSD')  <- getMeanSD(mask)
        
        if(cleanskips)
            mask <- cleanskips(mask)
    }

    mask
}

read.linearmask <- function (file = NULL, data = NULL, spacing = 10, spacingfactor = 1.5,
                             graph = TRUE, cleanskips = TRUE, ...)
{
    if (is.null(data) & !is.null(file)) {
        if (tools::file_ext(file) == 'shp')
            data <- readShapeSpatial(fn = tools::file_path_sans_ext(file))
        else
            data <- read.table (file, ...)
    }
    else if (is.null(data))
       stop("require one of 'file' or 'data'")

    isSLDF <- is(data, "SpatialLinesDataFrame")
    if (!isSLDF) {
        if (length(dim(data))!=2)
            stop ("require SpatialLinesDataFrame, dataframe or matrix",
                  " for 'data' input to read.linearmask")
        coln <- colnames(data)
        ixy <- match(c('x', 'y'), coln)
        if (any(is.na(ixy))) ixy <- 1:2
        mask <- as.data.frame(data[,ixy])
        names(mask) <- c('x', 'y')
        if (any(!apply(mask, 2, is.numeric)))
        stop ("non-numeric x or y coordinates")
        if (any(is.na(unlist(mask))))
            stop ("missing value(s) in x or y coordinates")
        if ('LineID' %in% coln)
            mask <- cbind(data[,'LineID'], mask)
        SLDF <- make.sldf(mask)
    }
    else {
        SLDF <- data
    }

    make.linearmask(SLDF, spacing, spacingfactor, graph, cleanskips)
}
