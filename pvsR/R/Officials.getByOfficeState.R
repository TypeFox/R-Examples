##' Get a list of officials according to office
##' 
##' This function is a wrapper for the Officials.getByOfficeState() method of the PVS API Officials class which grabs a list of officials according to office. The function sends a request with this method to the PVS API for all office IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Officials.getByOfficeState(stateId="NA", officeId)
##' @param stateId a character string or list of character strings with the state ID(s) (default is "NA", for national) (see references for details)
##' @param officeId a character string or list of character strings with the office ID(s) (see references for details)
##' @return A data frame with a row for each official and columns with the following variables describing the official:\cr candidateList.candidate*.candidateId,\cr candidateList.candidate*.firstName,\cr candidateList.candidate*.nickName,\cr candidateList.candidate*.middleName,\cr candidateList.candidate*.lastName,\cr candidateList.candidate*.suffix,\cr candidateList.candidate*.title,\cr candidateList.candidate*.electionParties,\cr candidateList.candidate*.officeParties,\cr candidatelist.candidate*.officeStatus,\cr candidateList.candidate*.officeDistrictId,\cr candidateList.candidate*.officeDistrictName,\cr candidateList.candidate*.officeTypeId,\cr candidateList.candidate*.officeId,\cr candidateList.candidate*.officeName,\cr candidateList.candidate*.officeStateId.
##' @references http://api.votesmart.org/docs/Officials.html\cr
##' Use State.getStateIDs() to get a list of state IDs.\cr
##' See http://api.votesmart.org/docs/semi-static.html or use Office.getOfficesByType(), Office.getOfficesByLevel(), Office.getOfficesByTypeLevel() or Office.getOfficesByBranchLevel() to get a list of office ID(s).
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of officials by state and office ID
##' \dontrun{officials <- Officials.getByOfficeState(,as.list(60:69))}
##' \dontrun{officials}
##' @export

Officials.getByOfficeState <-
function (stateId="NA", officeId) {
  
  
  # internal function
  Officials.getByOfficeState.basic <- function (.stateId, .officeId) {
    
    request <-  "Officials.getByOfficeState?"
    inputs  <-  paste("&stateId=",.stateId,"&officeId=",.officeId,sep="")
    output  <-  pvsRequest4(request,inputs)
    output$stateId <- .stateId
    output$officeId <- .officeId
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(stateId, FUN= function (y) {
    lapply(officeId, FUN= function (s) {
      Officials.getByOfficeState.basic(.stateId=y, .officeId=s)
    }
           )
  }
                        )
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
