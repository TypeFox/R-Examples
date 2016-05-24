#' Caluclate maximum lake length
#' 
#' Maximum lake length is defined as the longest open water distance of a lake.
#' This function determines the maximum lake length of lake by comparing the 
#' lenghts of a user defined number of lines.  The user specifies a number of 
#' points to distribute along the lake shoreline and the point to point line 
#' lengths are checked for multiple intersections (i.e. line not just open 
#' water),starting with the longest line first.  Function is currently very 
#' sensitive to \code{pointDens} and slows down considerably for larger values.
#' Small values of point dens are acceptable for lakes withtout complex 
#' shorelines.  Care needs to be taken in selecting an approriate 
#' \code{pointDens} and multiple values should be checked to ensure stability
#' in the estimates.  
#' 
#' @param inLakeMorpho An object of \code{\link{lakeMorphoClass}}.  Output of the 
#'        \code{\link{lakeSurroundTopo}} function would be appropriate as input
#' @param pointDens Number of points to place equidistant along shoreline. The 
#'        maximum point to point distance that does not also intersect the 
#'        shoreline is used.  To total of n*(n-1)/2 comparisons is possible, but
#'        in practice is usually significant less.
#' @param addLine Boolean to determine if the selected max length line should be 
#'        added to the inLakeMorpho object.  Defaults to True
#' 
#' @export
#' @return This returns a numeric value indicating the length of the longest 
#'         line possible in the lake. Units are the same as the input data.
#' 
#' @references Florida LAKEWATCH (2001). A Beginner's guide to water management
#'             - Lake Morphometry (2nd ed.). Gainesville: Florida LAKEWATCH, 
#'             Department of Fisheries and Aquatic Sciences.
#'             \href{http://edis.ifas.ufl.edu/pdffiles/FA/FA08100.pdf}{Link}
#' @import sp rgeos
#' @examples
#' data(lakes)
#' lakeMaxLength(inputLM,50)

lakeMaxLength <- function(inLakeMorpho, pointDens, addLine = T) {
    if (class(inLakeMorpho) != "lakeMorpho") {
        return(warning("Input data is not of class 'lakeMorpho'.  Run lakeSurround Topo first."))
    }
    result <- NA
    lakeShorePoints <- spsample(as(inLakeMorpho$lake, "SpatialLines"), pointDens, "regular")@coords
    dm <- dist(lakeShorePoints)
    md <- nrow(lakeShorePoints)
    x0 <- lakeShorePoints[which(lower.tri(matrix(1, md, md)) == 1, arr.ind = T)[, 1], ][, 1][order(dm, decreasing = T)]  #[30:md]
    y0 <- lakeShorePoints[which(lower.tri(matrix(1, md, md)) == 1, arr.ind = T)[, 1], ][, 2][order(dm, decreasing = T)]  #[30:md]
    x1 <- lakeShorePoints[which(lower.tri(matrix(1, md, md)) == 1, arr.ind = T)[, 2], ][, 1][order(dm, decreasing = T)]  #[30:md]
    y1 <- lakeShorePoints[which(lower.tri(matrix(1, md, md)) == 1, arr.ind = T)[, 2], ][, 2][order(dm, decreasing = T)]  #[30:md]
    xydf <- data.frame(x0, x1, y0, y1)
    xylist <- split(xydf, rownames(xydf))
    myLines <- SpatialLines(lapply(xylist, function(x) Lines(list(Line(matrix(as.numeric(x), 2, 2))), row.names(x))), 
        proj4string = CRS(proj4string(inLakeMorpho$lake)))
    myInd <- gWithin(myLines, inLakeMorpho$lake, byid = T)
    if (sum(myInd) == 0) {
        return(NA)
    }
    myLine <- myLines[myInd][gLength(myLines[myInd], byid = T) == max(gLength(myLines[myInd], byid = T))]
    result <- gLength(myLine)
    
    if (addLine) {
        myName <- paste(substitute(inLakeMorpho))
        inLakeMorpho$maxLengthLine <- NULL
        inLakeMorpho <- c(inLakeMorpho, maxLengthLine = myLine)
        class(inLakeMorpho) <- "lakeMorpho"
        assign(myName, inLakeMorpho, envir = parent.frame())
    }
    return(result)
} 
