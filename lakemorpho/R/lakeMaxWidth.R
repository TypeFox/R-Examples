#' Function to find line representing maximum lake Width
#' 
#' Maximum lake width is defined as the maximum in lake distance that is 
#' perpendicular to the maximum lake length.  This function calculates the 
#' equation of the perpendicular line and repeats that line \code{pointDens}
#' number of times and returns the longest of those lines.
#'
#' @param inLakeMorpho An object of \code{\link{lakeMorphoClass}}.  Output of the 
#'        \code{\link{lakeSurroundTopo}} function would be appropriate as input
#' @param pointDens Number of points to place equidistant along the 
#'        \code{\link{lakeMaxLength}}. A line that crosses at that point and
#'        extends from shore to shore is calcuated.   
#' @param addLine Boolean to determine if the selected max length line should be 
#'        added to the inLakeMorpho object.  Defaults to True
#' @export
#' @return Returns a numeric value indicating the length of the longest
#'         line perpndicular to the maximum length line.
#' 
#' @references Florida LAKEWATCH (2001). A Beginner's guide to water management
#'             - Lake Morphometry (2nd ed.). Gainesville: Florida LAKEWATCH, 
#'             Department of Fisheries and Aquatic Sciences.
#'             \href{http://edis.ifas.ufl.edu/pdffiles/FA/FA08100.pdf}{Link}
#' @import sp rgeos
#' 
#' @examples
#' data(lakes)
#' lakeMaxWidth(inputLM,50)

lakeMaxWidth <- function(inLakeMorpho, pointDens, addLine = T) {
    myName <- paste(substitute(inLakeMorpho))
    if (class(inLakeMorpho) != "lakeMorpho") {
        return(warning("Input data is not of class 'lakeMorpho'.  Run lakeSurround Topo first."))
    }
    if (is.null(inLakeMorpho$maxLengthLine)) {
        lakeMaxLength(inLakeMorpho, pointDens)
    }
    linedata <- data.frame(spsample(inLakeMorpho$maxLengthLine, 30, "regular")@coords)
    x <- linedata[, 1]
    y <- linedata[, 2]
    xdiff <- abs(x[1] - x[2])
    ydiff <- abs(y[1] - y[2])
    lakeMinx <- bbox(inLakeMorpho$lake)[1, 1]
    lakeMaxx <- bbox(inLakeMorpho$lake)[1, 2]
    lakeMiny <- bbox(inLakeMorpho$lake)[2, 1]
    lakeMaxy <- bbox(inLakeMorpho$lake)[2, 2]
    while (!(max(x) > lakeMaxx && min(x) < lakeMinx) && !(max(y) > lakeMaxy && min(y) < lakeMiny)) {
        if (x[1] - x[2] >= 0) {
            x <- c(x, x[length(x)] - xdiff)
            x <- c(x[1] + xdiff, x)
        } else {
            x <- c(x, x[length(x)] + xdiff)
            x <- c(x[1] - xdiff, x)
        }
        if (y[1] - y[2] >= 0) {
            y <- c(y, y[length(y)] - ydiff)
            y <- c(y[1] + ydiff, y)
        } else {
            y <- c(y, y[length(y)] + ydiff)
            y <- c(y[1] - ydiff, y)
        }
    }
    longline <- matrix(c(min(x), max(x), min(y), max(y)), 2, 2)
    longline <- SpatialLines(list(Lines(list(Line(longline)), "1")), proj4string = CRS(proj4string(inLakeMorpho$lake)))
    longlinedata <- data.frame(spsample(longline, pointDens, "regular")@coords)
    if (round(var(y), 1) == 0) {
        xmin <- longlinedata[, 1]
        xmax <- longlinedata[, 1]
        ypred_min <- bbox(inLakeMorpho$lake)[2, 1]
        ypred_max <- bbox(inLakeMorpho$lake)[2, 2]
    } else if (round(var(x), 1) == 0) {
        xmin <- bbox(inLakeMorpho$lake)[1, 1]
        xmax <- bbox(inLakeMorpho$lake)[1, 2]
        ypred_min <- longlinedata[, 2]
        ypred_max <- longlinedata[, 2]
    } else {
        mylm <- lm(longlinedata[, 2] ~ longlinedata[, 1])
        mylm2_slope <- (1/mylm$coefficients[2]) * -1
        mylm2_int <- (mylm2_slope * -longlinedata[, 1]) - (-longlinedata[, 2])
        xmin <- bbox(inLakeMorpho$lake)[1, 1]
        xmax <- bbox(inLakeMorpho$lake)[1, 2]
        ypred_min <- (mylm2_slope * xmin) + mylm2_int
        ypred_max <- (mylm2_slope * xmax) + mylm2_int
    }
    mydf <- data.frame(xmin, xmax, ypred_min, ypred_max, ID = 1:length(longlinedata[, 2]))
    createSL <- function(x) {
        mat <- matrix(as.numeric(x[1:4]), 2, 2)
        id <- as.numeric(x[5])
        return(Lines(list(Line(mat)), id))
    }
    mylinelist <- apply(mydf, 1, createSL)
    mylines <- SpatialLines(mylinelist, proj4string = CRS(proj4string(inLakeMorpho$lake)))
    myInter <- gIntersection(mylines[gCrosses(mylines, inLakeMorpho$lake, byid = T), ], inLakeMorpho$lake, 
        byid = T)
    lineInter <- unlist(lapply(myInter@lines, function(x) slot(x, "Lines")))
    myInter2 <- list()
    for (i in 1:length(lineInter)) {
        myInter2 <- c(myInter2, Lines(lineInter[i], i))
    }
    myInter2Lines <- SpatialLines(myInter2, proj4string = CRS(proj4string(inLakeMorpho$lake)))
    maxWidthLine <- myInter2Lines[gLength(myInter2Lines, byid = T) == max(gLength(myInter2Lines, byid = T)), 
        ]
    if (addLine) {
        
        inLakeMorpho$maxWidthLine <- NULL
        inLakeMorpho <- c(inLakeMorpho, maxWidthLine = maxWidthLine)
        class(inLakeMorpho) <- "lakeMorpho"
        assign(myName, inLakeMorpho, envir = parent.frame())
    }
    return(gLength(maxWidthLine))
} 
