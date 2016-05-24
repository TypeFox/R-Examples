##' Get a list of special interest groups according to rating category and state.
##' 
##' This function is a wrapper for the Rating.getSigList() method of the PVS API Rating class which dumps special interest groups (SIGs) according to category and state. The function sends a request with this method to the PVS API for all category and state IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Rating.getSigList(stateId="NA", categoryId)
##' @param stateId (optional) a character string or list of character strings with the state ID(s) (default: "NA", for national) (see references for details)
##' @param categoryId a character string or list of character strings with the category ID(s) (see references for details)
##' @return A data frame with a row for each special interest group and columns with the following variables describing the special interest group:\cr sigs.sig*.sigId,\cr sigs.sig*.parentId,\cr sigs.sig*.name.
##' @references http://api.votesmart.org/docs/Rating.html\cr
##' Use State.getStateIDs() to get a list of state IDs.\cr
##' Use Rating.getCategories() or Rating.getCandidateRating() to get category IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of special interest groups for certain categories and all states
##' \dontrun{sig <- Rating.getSigList(categoryId=list(2,4,5,7))}
##' \dontrun{sig}
##' @export



Rating.getSigList <-
function (stateId="NA", categoryId) {
  
  
  # internal function
  Rating.getSigList.basic <- function (.stateId, .categoryId) {
    
    request <-  "Rating.getSigList?"
    inputs  <-  paste("&stateId=",.stateId,"&categoryId=",.categoryId,sep="")
    output  <-  pvsRequest4(request,inputs)
    output$stateId <- .stateId
    output$categoryId <- .categoryId
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(stateId, FUN= function (y) {
    lapply(categoryId, FUN= function (s) {
      Rating.getSigList.basic(.stateId=y, .categoryId=s)
    }
           )
  }
                        )
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
