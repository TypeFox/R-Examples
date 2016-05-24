#' Function to calculate fetch along an input bearing
#' 
#' The function calculates the maximum in lake distance of a line along an input
#' bearing.
#' 
#' @param inLakeMorpho An object of \code{\link{lakeMorphoClass}}.  Output of the 
#'        \code{\link{lakeSurroundTopo}} function would be appropriate as input
#' @param bearing Character that indicates the bearing of the desired fetch 
#' @param addLine Boolean to determine if the selected max length line should be 
#'        added to the inLakeMorpho object.  Defaults to True.  Note that the 
#'        line is returned in the same projection as the input data.  
#' @export
#' @return Returns a numeric value indicating the length of the longest 
#'         line in the lake along the input bearing. Units are the same as the 
#'         input data.
#' 
#' @references Florida LAKEWATCH (2001). A Beginner's guide to water management
#'             - Lake Morphometry (2nd ed.). Gainesville: Florida LAKEWATCH, 
#'             Department of Fisheries and Aquatic Sciences.
#'             \href{http://edis.ifas.ufl.edu/pdffiles/FA/FA08100.pdf}{Link}
#' 
#' @import sp geosphere rgeos rgdal maptools
#' 
#' @export
#' 
#' @examples
#' data(lakes)
#' lakeFetch(inputLM,45)

lakeFetch <- function(inLakeMorpho, bearing, addLine = T) {
    inputName <- paste(substitute(inLakeMorpho))
    if (class(inLakeMorpho) != "lakeMorpho") {
        return(warning("Input data is not of class 'lakeMorpho'.  Run lakeSurround Topo first."))
    }
    result <- NA
    # convert to dd
    lakedd <- spTransform(inLakeMorpho$lake, CRSobj = CRS("+proj=longlat +datum=WGS84"))
    # get min/max distance: converts original extent to square.  ensures full coverage of possible lines
    origMinMin <- SpatialPoints(matrix(bbox(lakedd)[, 1], 1, 2), proj4string = CRS("+proj=longlat +datum=WGS84"))
    origMaxMax <- SpatialPoints(matrix(bbox(lakedd)[, 2], 1, 2), proj4string = CRS("+proj=longlat +datum=WGS84"))
    origMinMax <- SpatialPoints(matrix(c(bbox(lakedd)[1, 1], bbox(lakedd)[2, 2]), 1, 2), proj4string = CRS("+proj=longlat +datum=WGS84"))
    origMaxMin <- SpatialPoints(matrix(c(bbox(lakedd)[1, 2], bbox(lakedd)[2, 1]), 1, 2), proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    # Get distances for each side of bounding box
    l1 <- distCosine(origMinMin, origMaxMin)
    l2 <- distCosine(origMinMin, origMinMax)
    # get new points to make the extent square
    if (l1 > l2) {
        minPt <- SpatialPoints(destPoint(origMinMin, 180, (l1 - l2)/2), proj4string = CRS("+proj=longlat +datum=WGS84"))
        maxPt <- SpatialPoints(destPoint(origMaxMax, 0, (l1 - l2)/2), proj4string = CRS("+proj=longlat +datum=WGS84"))
    } else {
        minPt <- SpatialPoints(destPoint(origMinMin, 270, (l2 - l1)/2), proj4string = CRS("+proj=longlat +datum=WGS84"))
        maxPt <- SpatialPoints(destPoint(origMaxMax, 90, (l2 - l1)/2), proj4string = CRS("+proj=longlat +datum=WGS84"))
    }
    maxDist <- distCosine(minPt, maxPt)
    
    # Convert bearing to half
    if (bearing > 180) {
        bearing <- bearing - 180
    }
    # calc perpendicular bearings
    if (bearing < 90) {
        perpbear1 <- bearing + 90
        perpbear2 <- perpbear1 + 180
    } else {
        perpbear1 <- bearing - 90
        perpbear2 <- perpbear1 + 180
    }
    
    # Build list of center points for perpbear1
    centPts <- list()
    centPts[[1]] <- coordinates(lakedd)
    colnames(centPts[[1]]) <- c("lon", "lat")
    centPts[[2]] <- destPoint(centPts[[1]], perpbear1, max(res(inLakeMorpho$lakeDistance)) * 3)
    i <- length(centPts)
    while (centPts[[i]][, 1] < coordinates(maxPt)[, 1] & centPts[[i]][, 1] > coordinates(minPt)[, 1] & centPts[[i]][, 
        2] < coordinates(maxPt)[, 2] & centPts[[i]][, 2] > coordinates(minPt)[, 2]) {
        i <- length(centPts) + 1
        centPts[[i]] <- destPoint(centPts[[i - 1]], perpbear1, round(max(res(inLakeMorpho$lakeDistance)) * 
            3))
    }
    # Build list of center points for perpbear2
    i <- length(centPts) + 1
    centPts[[i]] <- destPoint(centPts[[1]], perpbear2, max(res(inLakeMorpho$lakeDistance)) * 3)
    while (centPts[[i]][, 1] < coordinates(maxPt)[, 1] & centPts[[i]][, 1] > coordinates(minPt)[, 1] & centPts[[i]][, 
        2] < coordinates(maxPt)[, 2] & centPts[[i]][, 2] > coordinates(minPt)[, 2]) {
        i <- length(centPts) + 1
        centPts[[i]] <- destPoint(centPts[[i - 1]], perpbear2, round(max(res(inLakeMorpho$lakeDistance)) * 
            3))
    }
    # calc point for centroid, max distance, bearing + 180 (if bearing is less that 180) or - 180 (if bearing
    # is more than 180)
    allLines <- list()
    for (i in 1:length(centPts)) {
        allLines[[i]] <- Lines(list(Line(rbind(destPoint(centPts[[i]], bearing, maxDist), destPoint(centPts[[i]], 
            bearing + 180, maxDist)))), as.character(i))
    }
    allLinesSL <- SpatialLines(allLines, proj4string = CRS("+proj=longlat +datum=WGS84"))
    
    # clip out lines that are inside lake
    lakeLinesSL <- gIntersection(lakedd, allLinesSL, byid = TRUE)
    
    myInd <- gWithin(lakeLinesSL, lakedd, byid = T)
    lakeLinesSL_proj <- spTransform(lakeLinesSL, CRSobj = CRS(proj4string(inLakeMorpho$lake)))
    # Loop through each item in lakeLinesSL_Proj and create a SpatialLines object for each segment.
    lakeLinesList_proj <- list()
    for (i in 1:length(lakeLinesSL_proj)) {
        xlines <- slot(lakeLinesSL_proj[i], "lines")
        xLines <- slot(xlines[[1]], "Lines")
        for (j in 1:length(xLines)) {
            lakeLinesList_proj[length(lakeLinesList_proj) + 1] <- Lines(xLines[j], as.character(length(lakeLinesList_proj) + 
                1))
        }
    }
    
    # Determine the longest
    lakeLinesSL_proj <- SpatialLines(lakeLinesList_proj, proj4string = CRS(proj4string(inLakeMorpho$lake)))
    result <- max(gLength(lakeLinesSL_proj, byid = TRUE), na.rm = T)
    myLine <- lakeLinesSL_proj[gLength(lakeLinesSL_proj, byid = T) == result, ]
    
    # line added to input lakemorpho
    if (addLine) {
        myName <- paste("maxFetchLine_", bearing, sep = "")
        inLakeMorpho[[substitute(myName)]] <- NULL
        inLakeMorpho[[substitute(myName)]] <- myLine
        class(inLakeMorpho) <- "lakeMorpho"
        assign(inputName, inLakeMorpho, envir = parent.frame())
    }
    return(result)
} 
