##' Get a list of officials according to state representation
##' 
##' This function is a wrapper for the Officials.getStatewide() method of the PVS API Officials class which grabs a list of officials according to the state they are representing. The function sends a request with this method to the PVS API for all state IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Officials.getStatewide(stateId="NA")
##' @param stateId (optional) a character string or list of character strings with the state ID(s) (default: "NA", for national) (see references for details)
##' @return A data frame with a row for each official and columns with the following variables describing the official:\cr candidateList.candidate*.candidateId,\cr candidateList.candidate*.firstName,\cr candidateList.candidate*.nickName,\cr candidateList.candidate*.middleName,\cr candidateList.candidate*.lastName,\cr candidateList.candidate*.suffix,\cr candidateList.candidate*.title,\cr candidateList.candidate*.electionParties,\cr candidateList.candidate*.officeParties,\cr candidatelist.candidate*.officeStatus,\cr candidateList.candidate*.officeDistrictId,\cr candidateList.candidate*.officeDistrictName,\cr candidateList.candidate*.officeTypeId,\cr candidateList.candidate*.officeId,\cr candidateList.candidate*.officeName,\cr candidateList.candidate*.officeStateId.
##' @references http://api.votesmart.org/docs/Officials.html\cr
##' Use State.getStateIDs() to get a list of state IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get all officials of a certain state
##' \dontrun{officials <- Officials.getStatewide("FL")}
##' \dontrun{head(officials)}

##' @export



Officials.getStatewide <-
function (stateId="NA") {
  
  
  # internal function
  Officials.getStatewide.basic <- function (.stateId) {
    
    request <-  "Officials.getStatewide?"
    inputs  <-  paste("&stateId=",.stateId,sep="")
    output  <-  pvsRequest4(request,inputs)
    output$stateId <- .stateId
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(stateId, FUN= function (s) {
    Officials.getStatewide.basic(.stateId=s)
  }
                        )
  
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
