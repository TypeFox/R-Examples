##' Get detailed information about a special interest group
##' 
##' This function is a wrapper for the Rating.getSig() method of the PVS API Rating class which dumps detailed information about special interest groups (SIGs). The function sends a request with this method to the PVS API for all SIG IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Rating.getSig(sigId)
##' @param sigId a character string or list of character strings with the special interest group's ID(s) (see references for details)
##' @return A data frame with a row for each special interest group and columns with the following variables describing the special interest group:\cr sig.sigId,\cr sig.parentId,\cr sig.stateId,\cr sig.name,\cr sig.description,\cr sig.address,\cr sig.city,\cr sig.state,\cr sig.zip,\cr sig.phone1,\cr sig.phone2,\cr sig.fax,\cr sig.email,\cr sig.url,\cr sig.contactName.
##' @references http://api.votesmart.org/docs/Rating.html\cr
##' Use Rating.getSigList() to get a list of special interest group's IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get information about certain special interest groups
##' \dontrun{info <- Rating.getSig(list(1016,1120))}
##' \dontrun{info}

##' @export










Rating.getSig <-
function (sigId) {
  
  
  # internal function
  Rating.getSig.basic <- function (.sigId) {
    
    request <-  "Rating.getSig?"
    inputs  <-  paste("&sigId=",.sigId,sep="")
    output  <-  pvsRequest4.1(request,inputs)
    output$sigId <- .sigId
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(sigId, FUN= function (s) {
    Rating.getSig.basic(.sigId=s)
  }
                        )
  
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
