getXYmat <-
function(DF, xyTopLeft=TRUE, relPOA=TRUE) {
    if(!is.data.frame(DF)) { stop("DF must be a data frame") }

    ## convert DF names to lower case
    DF <- setNames(DF, tolower(names(DF)))

    ## make sure DF has the required variable names
    ## coordinates need to be called either X, Y order Point.X, Point.Y
    dfNames  <- names(DF)                # what variables are present
    needsXY1 <- c("point.x", "point.y")  # coordinates must have this name
    needsXY2 <- c("x", "y")              # or this
    wantsAIM <- c("aim.x", "aim.y")      # useful
    hasXY1   <- needsXY1 %in% dfNames
    hasXY2   <- needsXY2 %in% dfNames
    hasAIM   <- wantsAIM %in% dfNames    # useful ones we have

    if(!xor(all(hasXY1), all(hasXY2))) {
       stop("Coordinates must be named either X, Y or Point.X, Point.Y")
    }

    if(("z" %in% dfNames) && ("point.z" %in% dfNames)) {
        stop("Coordinates must be named either Z or Point.Z")
    }

    ## analysis should be relative to POA, but POA is missing
    if(!all(hasAIM) && relPOA) {
        warning(c("The data frame is missing variable(s)\n",
                  paste(wantsAIM[!hasAIM], collapse=" "), "\n",
                  "Point of Aim is assumed to be in (0,0)"))
        relPOA <- FALSE
    }

    if(!relPOA) {                        # coords not relative to POA
        DF$aim.x <- 0                    # -> set POA to (0,0)
        DF$aim.y <- 0

        if(("z" %in% dfNames) || ("point.z" %in% dfNames)) {
            DF$aim.z <- 0
        }
    }
    
    ## if names are X, Y, Z rename to Point.X, Point.Y, Point.Z
    if(all(hasXY2)) {
        dfNames <- names(DF)
        dfNames[dfNames %in% "x"] <- "point.x"
        dfNames[dfNames %in% "y"] <- "point.y"
        dfNames[dfNames %in% "z"] <- "point.z"
        warning("Variables X, Y were renamed to Point.X, Point.Y")
        names(DF) <- dfNames
    }
    
    ## coords relative to point of aim
    ## y-coords exported from OnTarget: (0,0) is top-left
    x <- DF$point.x - DF$aim.x           # x-coords
    y <- if(xyTopLeft) {
        -(DF$point.y - DF$aim.y)         # y-coords
    } else {
          DF$point.y - DF$aim.y
    }

    z <- if(("z" %in% dfNames) || ("point.z" %in% dfNames)) {
        DF$point.z - DF$aim.z            # z-coords
    } else {
        NULL
    }

    return(cbind(x, y, z))               # new (x,y)-coords as matrix
}
