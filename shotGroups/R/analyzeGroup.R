analyzeGroup <-
function(DF, xyTopLeft=TRUE, conversion="m2cm", bandW=0.5,
         CEPtype="CorrNormal", bootCI=c("basic", "bca")) {
    if(!is.data.frame(DF)) { stop("DF must be a data frame") }

    ## convert DF names to lower case
    DF <- setNames(DF, tolower(names(DF)))

    #####-----------------------------------------------------------------------
    ## make sure DF has the required variable names
    varNames <- names(DF)                # what variables are present
    needsXY1 <- c("point.x", "point.y")  # coordinates must have this name
    needsXY2 <- c("x", "y")              # or this
    wantsDst <- "distance"               # useful
    wantsAIM <- c("aim.x", "aim.y")      # useful
    hasXY1   <- needsXY1 %in% varNames
    hasXY2   <- needsXY2 %in% varNames
    hasDst   <- wantsDst %in% varNames   # useful ones we have
    hasAIM   <- wantsAIM %in% varNames   # useful ones we have

    if(!xor(all(hasXY1), all(hasXY2))) {
        stop("Coordinates must be named either X, Y or Point.X, Point.Y")
    }

    if(!all(hasDst)) {
        warning(c("The data frame is missing variable\n",
                  paste(wantsDst[!hasDst], collapse=" "), "\n",
                  "Distance is assumed to be 100"))
        DF$distance <- 100
    }

    if(!all(hasAIM)) {
        warning(c("The data frame is missing variable(s)\n",
            paste(wantsAIM[!hasAIM], collapse=" "), "\n",
            "Point of Aim is assumed to be in (0,0)"))
        DF$aim.x <- 0
        DF$aim.y <- 0
    }

    #####-----------------------------------------------------------------------
    ## prepare data: get (x,y)-coords relative to point of aim as matrix
    xy <- getXYmat(DF, xyTopLeft=xyTopLeft)

    #####-----------------------------------------------------------------------
    ## assess shape, location and spread
    if(length(unique(DF$distance)) > 1L) {
        warning("Distance to target is not homogeneous - using the average distance")
    }

    dstTrgt  <- mean(DF$distance)       # distance to target
    shape    <- groupShape(xy, plots=TRUE, bandW=bandW, outlier="mcd",
                           dstTarget=dstTrgt, conversion=conversion)
    location <- groupLocation(xy, plots=FALSE, bootCI=bootCI, level=0.95,
                              dstTarget=dstTrgt, conversion=conversion)
    spread   <- groupSpread(xy, plots=TRUE, CEPlevel=0.5, CIlevel=0.95,
                            CEPtype=CEPtype, bootCI=bootCI,
                            dstTarget=dstTrgt, conversion=conversion)

    #####-----------------------------------------------------------------------
    ## return all the collected numerical results and tests
    return(c(shape, location, spread))
}
