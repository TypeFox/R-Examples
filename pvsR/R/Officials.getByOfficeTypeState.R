##' Get a list of officials according to office type and state
##' 
##' This function is a wrapper for the Officials.getByOfficeTypeState() method of the PVS API Officials class which grabs a list of officials according to the office type and state they represent. The function sends a request with this method to the PVS API for all state and office type IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Officials.getByOfficeTypeState(stateId="NA", officeTypeId)
##' @param stateId (optional) a character string or list of character strings with the state ID(s) (default: "NA", for national) (see references for details)
##' @param officeTypeId a character string or list of character strings with the office type ID(s) (see references for details)
##' @return A data frame with a row for each official and columns with the following variables describing the official:\cr candidateList.candidate*.candidateId,\cr candidateList.candidate*.firstName,\cr candidateList.candidate*.nickName,\cr candidateList.candidate*.middleName,\cr candidateList.candidate*.lastName,\cr candidateList.candidate*.suffix,\cr candidateList.candidate*.title,\cr candidateList.candidate*.electionParties,\cr candidateList.candidate*.officeParties,\cr candidatelist.candidate*.officeStatus,\cr candidateList.candidate*.officeDistrictId,\cr candidateList.candidate*.officeDistrictName,\cr candidateList.candidate*.officeTypeId,\cr candidateList.candidate*.officeId,\cr candidateList.candidate*.officeName,\cr candidateList.candidate*.officeStateId.
##' @references http://api.votesmart.org/docs/Officials.html\cr
##' Use State.getStateIDs() to get a list of state IDs.\cr
##' See http://api.votesmart.org/docs/semi-static.html or use Office.getTypes or Office.getOfficesByLevel to get a list of office types ID(s). 
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' # Note that some officeTypeIds are only available on the state level or national level 
##' # (e.g. "L" for State Legislature only if stateId is specified!)
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of officials by state and office type 
##' \dontrun{CAlegislators <- Officials.getByOfficeTypeState(officeTypeId="L", stateId="CA")}
##' \dontrun{head(CAlegislators)}
##' \dontrun{suprcourt <- Officials.getByOfficeTypeState(officeTypeId="J")}
##' \dontrun{head(suprcourt)}
##' @export



Officials.getByOfficeTypeState <-
function (stateId="NA", officeTypeId) {
  
  
  # internal function
  Officials.getByOfficeTypeState.basic <- function (.stateId, .officeTypeId) {
    
    request <-  "Officials.getByOfficeTypeState?"
    inputs  <-  paste("&stateId=",.stateId,"&officeTypeId=",.officeTypeId,sep="")
    output  <-  pvsRequest4(request,inputs)
    output$stateId <- .stateId
    output$officeTypeId <- .officeTypeId
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(stateId, FUN= function (y) {
    lapply(officeTypeId, FUN= function (s) {
      Officials.getByOfficeTypeState.basic(.stateId=y, .officeTypeId=s)
    }
           )
  }
                        )
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
}
