#' Return lake surface area
#' 
#' This function simply retruns the area of the lake SpatialPolygons that is part
#' of the \code{\link{lakeMorphoClass}} class.
#' 
#' @param inLakeMorpho an object of \code{\link{lakeMorphoClass}}.  Output of the 
#'        \code{\link{lakeSurroundTopo}} function would be appropriate as input
#'          
#' @export
#' @return Returns a numeric value for the total surface area of the lake
#' @references Florida LAKEWATCH (2001). A Beginner's guide to water management
#'             - Lake Morphometry (2nd ed.). Gainesville: Florida LAKEWATCH, 
#'             Department of Fisheries and Aquatic Sciences.
#'             \href{http://edis.ifas.ufl.edu/pdffiles/FA/FA08100.pdf}{Link}
#' @import rgeos
#' @examples
#' data(lakes)
#' lakeSurfaceArea(inputLM)

lakeSurfaceArea <- function(inLakeMorpho) {
    if (class(inLakeMorpho) != "lakeMorpho") {
        return(warning("Input data is not of class 'lakeMorpho'.  Run lakeSurround Topo first."))
    }
    return(gArea(inLakeMorpho$lake))
} 
