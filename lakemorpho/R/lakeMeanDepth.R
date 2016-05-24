#' Function to return lake Mean Depth
#' 
#' Calculates average depth of lake as a mean of lake volume divided by
#' lake surface area
#' 
#' @param inLakeMorpho An object of \code{\link{lakeMorphoClass}}.  Output of the 
#'        \code{\link{lakeSurroundTopo}} function would be appropriate as input
#' @export      
#' @return Returns a numerica value for the mean depth of the lake 
#' @references Florida LAKEWATCH (2001). A Beginner's guide to water management
#'             - Lake Morphometry (2nd ed.). Gainesville: Florida LAKEWATCH, 
#'             Department of Fisheries and Aquatic Sciences.
#'             \href{http://edis.ifas.ufl.edu/pdffiles/FA/FA08100.pdf}{Link}
#' 
#' @examples
#' data(lakes)
#' lakeMeanDepth(inputLM)


lakeMeanDepth <- function(inLakeMorpho) {
    if (class(inLakeMorpho) != "lakeMorpho") {
        return(warning("Input data is not of class 'lakeMorpho'.  Run lakeSurround Topo first."))
    }
    return(lakeVolume(inLakeMorpho)/lakeSurfaceArea(inLakeMorpho))
} 
