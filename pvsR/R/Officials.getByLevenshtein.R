##' Get a list of officials according to an approximate last name match
##' 
##' This function is a wrapper for the Officials.getByLevenshtein() method of the PVS API Officials class which grabs a list of officials according to an approximate last name match. The function sends a request with this method to the PVS API for all last names given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Officials.getByLevenshtein(lastName)
##' @param lastName a character string or list of character strings with the last name(s) (see references for details)
##' @return A data frame with a row for each official and columns with the following variables describing the official:\cr candidateList.candidate*.candidateId,\cr candidateList.candidate*.firstName,\cr candidateList.candidate*.nickName,\cr candidateList.candidate*.middleName,\cr candidateList.candidate*.lastName,\cr candidateList.candidate*.suffix,\cr candidateList.candidate*.title,\cr candidateList.candidate*.electionParties,\cr candidateList.candidate*.officeParties,\cr candidatelist.candidate*.officeStatus,\cr candidateList.candidate*.officeDistrictId,\cr candidateList.candidate*.officeDistrictName,\cr candidateList.candidate*.officeTypeId,\cr candidateList.candidate*.officeId,\cr candidateList.candidate*.officeName,\cr candidateList.candidate*.officeStateId.
##' @references http://api.votesmart.org/docs/Officials.html
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of officials with similar last names
##' \dontrun{names <- Officials.getByLevenshtein(list("Miller","Fine"))}
##' \dontrun{head(names)}

##' @export

Officials.getByLevenshtein <-
function (lastName) {
  
  
  # internal function
  Officials.getByLevenshtein.basic <- function (.lastName) {
    
    request <-  "Officials.getByLevenshtein?"
    inputs  <-  paste("&lastName=",.lastName,sep="")
    output  <-  pvsRequest4(request,inputs)
    output$lastName.input <- .lastName
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(lastName, FUN= function (y) {
    
    Officials.getByLevenshtein.basic(.lastName=y)
    
  }
                        )
  
  
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
