#' Calculate all Lake Morphometry Metrics
#' 
#' This function is a wrapper function for all of the metrics. It calculates
#' each metric for an input lakeMorphoClass.  This returns a list of all 
#' metrics
#' 
#' @param inLakeMorpho an object of \code{\link{lakeMorphoClass}}.  Output of the 
#'        \code{\link{lakeSurroundTopo}} function would be appropriate as input
#' @param bearing Character that indicates the bearing of the desired fetch
#' @param pointDens Number of points to place equidistant along shoreline for 
#'        \code{\link{lakeMaxLength}} or density of lines to test for 
#'        \code{\link{lakeMaxWidth}} and \code{\link{lakeFetch}}.
#' @param correctFactor Value used to correct the predicted maximum lake depth.  
#'        Defaults to 1. Corrections are simply accomplished by multiplying 
#'        estimated max depth by correction factor. Correction factors can be 
#'        determined empirically by regressing the predicted depth against a 
#'        known maximum depth while forcing the intercept through zero.  The 
#'        slope of the line would then be used as the correction 
#'        factor(Hollister et. al, 2011). 
#'        
#'         
#' @references Florida LAKEWATCH (2001). A Beginner's guide to water management
#'             - Lake Morphometry (2nd ed.). Gainesville: Florida LAKEWATCH, 
#'             Department of Fisheries and Aquatic Sciences.
#'             \href{http://edis.ifas.ufl.edu/pdffiles/FA/FA08100.pdf}{Link}
#' @references Hollister, J. W., W.B. Milstead (2010). Using GIS to Estimate 
#'             Lake Volume from Limited Data. Lake and Reservoir Management. 
#'             26(3)194-199.
#'             \href{http://dx.doi.org/10.1080/07438141.2010.504321}{Link}
#' @references Hollister, J. W., W.B. Milstead, M.A. Urrutia (2011). Predicting 
#'             Maximum Lake Depth from Surrounding Topography. PLoS ONE 6(9).
#'             \href{http://dx.doi.org/10.1371/journal.pone.0025764}{link} 
#' @export
#' @return Returns a list with all lake metrics calculated for a given input
#'         lakemorpho object
#' 
#' @examples
#' \dontrun{
#' data(lakes)
#' calcLakeMetrics(inputLM,45,250)}

calcLakeMetrics <- function(inLakeMorpho, bearing, pointDens, correctFactor = 1) {
    if (class(inLakeMorpho) != "lakeMorpho") {
        return(warning("Input data is not of class 'lakeMorpho'.  Run lakeSurround Topo first."))
    }
    allMet <- list(surfaceArea = lakeSurfaceArea(inLakeMorpho), shorelineLength = lakeShorelineLength(inLakeMorpho), 
        shorelinDevelopment = lakeShorelineDevelopment(inLakeMorpho), maxDepth = lakeMaxDepth(inLakeMorpho, 
            correctFactor), volume = lakeVolume(inLakeMorpho, correctFactor), meanDepth = lakeMeanDepth(inLakeMorpho), 
        maxLength = lakeMaxLength(inLakeMorpho, pointDens), maxWidth = lakeMaxWidth(inLakeMorpho, pointDens), 
        meanWidth = lakeMeanWidth(inLakeMorpho), fetch = lakeFetch(inLakeMorpho, bearing))
    return(allMet)
} 
