#' Get details for a given feature
#' 
#' Given the featureID (an integer), this function get details
#' of the taxon from the NBN Gateway. These include the label (name),
#' feature type, etc. If the feature is a \code{GridSquare} then
#' the information includes the bounding box coordinates.
#' 
#' @export
#' @param featureID A the featureID as a string
#' @return A list containing the JSON object returned by the NBN Gateway.
#' @author Stuart Ball, JNCC \email{stuart.ball@@jncc.gov.uk}
#' @seealso \code{\link{getOccurrences}}
#' @examples \dontrun{ 
#'  t <- getFeature("97479")
#'  t['label']  ## [1] "SN413499"
#' }
#' 
getFeature <- function(featureID=NULL) {
    runnbnurl(service="feature", feature=featureID)
}
