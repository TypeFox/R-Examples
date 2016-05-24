##' Get categories that contain released ratings according to state
##' 
##' This function is a wrapper for the Rating.getCategories() method of the PVS API Rating class which dumps categories that contain released ratings according to state. The function sends a request with this method to the PVS API for all state IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Rating.getCategories(stateId="NA")
##' @param stateId (optional) a character string or list of character strings with the stateId(s) (default: "NA", for national) (see references for details)
##' @return A data frame with a row for each rating and columns with the following variables describing the rating:\cr categories.category*.categoryId,\cr categories.category*.name.
##' @references http://api.votesmart.org/docs/Rating.html\cr
##' Use State.getStateIDs() to get a list of state IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get rating categories
##' \dontrun{rating <- Rating.getCategories()}
##' \dontrun{rating}

##' @export







Rating.getCategories <-
function (stateId="NA") {
  
  
  # internal function
  Rating.getCategories.basic <- function (.stateId) {
    
    request <-  "Rating.getCategories?"
    inputs  <-  paste("&stateId=",.stateId,sep="")
    output  <-  pvsRequest4(request,inputs)
    output$stateId <- .stateId
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(stateId, FUN= function (s) {
    Rating.getCategories.basic(.stateId=s)
  }
                        )
  
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
