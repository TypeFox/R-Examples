##' Get all candidate ratings from an evaluation by a special interest group
##' 
##' This function is a wrapper for the Rating.getRating() method of the PVS API Rating class which dumps all candidate ratings from a scorecard by a special interest group (SIG). The function sends a request with this method to the PVS API for all rating IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Rating.getRating(ratingId)
##' @param ratingId a character string or list of character strings with the rating ID(s) (see references for details)
##' @return A data frame with a row for each candidate and columns with the following variables describing the candidate:\cr candidateRating*.candidateId,\cr candidateRating*.rating.
##' @references http://api.votesmart.org/docs/Rating.html\cr
##' Use Rating.getSigRatings() or Rating.getCandidateRating() to get a list of rating IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get the candidate rating by a certain special interest group
##' \dontrun{scorecard <- Rating.getRating(77)}
##' \dontrun{scorecard}

##' @export



Rating.getRating <-
function (ratingId) {
  
  
  # internal function
  Rating.getRating.basic <- function (.ratingId) {
    
    request <-  "Rating.getRating?"
    inputs  <-  paste("&ratingId=",.ratingId,sep="")
    output  <-  pvsRequest4(request,inputs)
    output$ratingId <- .ratingId
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(ratingId, FUN= function (s) {
    Rating.getRating.basic(.ratingId=s)
  }
                        )
  
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
