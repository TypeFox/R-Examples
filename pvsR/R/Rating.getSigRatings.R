##' Get all ratings (scorecards) by a special interest group
##' 
##' This function is a wrapper for the Rating.getSigRatings() method of the PVS API Rating class which dumps all ratings (scorecards) by a special interest group. The function sends a request with this method to the PVS API for all candidate and SIG IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Rating.getSigRatings(sigId)
##' @param sigId a character string or list of character strings with the special interest group's ID(s) (see references for details)
##' @return A data frame with a row for each special interest group and columns with the following variables describing the rating:\cr sig.sigId,\cr sig.name,\cr rating*.ratingId,\cr rating*.timespan,\cr rating*.ratingName,\cr rating*.ratingText.
##' @references http://api.votesmart.org/docs/Rating.html\cr
##' Use Rating.getSigList() to get a list of special interest group's IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get ratings of certain special interest groups
##' \dontrun{rating <- Rating.getSigRatings(list(568,1120,1704))}
##' \dontrun{rating}

##' @export





Rating.getSigRatings <-
function (sigId) {
  
  
  # internal function
  Rating.getSigRatings.basic <- function (.sigId) {
    
    request <-  "Rating.getSigRatings?"
    inputs  <-  paste("&sigId=",.sigId,sep="")
    output  <-  pvsRequest7(request,inputs)
    output$sigId <- .sigId
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(sigId, FUN= function (s) {
    Rating.getSigRatings.basic(.sigId=s)
  }
                        )
  
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
