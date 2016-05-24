##' Get a list of officials according to the district they are running for
##' 
##' This function is a wrapper for the Officials.getByDistrict() method of the PVS API Officials class which grabs a list of officials according to the district they are running for. The function sends a request with this method to the PVS API for all district IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Officials.getByDistrict(districtId)
##' @param districtId a character string or list of character strings with the district ID(s) (see references for details)
##' @return A data frame with a row for each official and columns with the following variables describing the official:\cr candidateList.candidate*.candidateId,\cr candidateList.candidate*.firstName,\cr candidateList.candidate*.nickName,\cr candidateList.candidate*.middleName,\cr candidateList.candidate*.lastName,\cr candidateList.candidate*.suffix,\cr candidateList.candidate*.title,\cr candidateList.candidate*.electionParties,\cr candidateList.candidate*.electionstatus,\cr candidateList.candidate*.officeParties,\cr candidatelist.candidate*.officeStatus,\cr candidateList.candidate*.officeDistrictId,\cr candidateList.candidate*.officeDistrictName,\cr candidateList.candidate*.officeTypeId,\cr candidateList.candidate*.officeId,\cr candidateList.candidate*.officeName,\cr candidateList.candidate*.officeStateId.
##' @references http://api.votesmart.org/docs/Officials.html\cr
##' Use District.getByOfficeState() or District.getByZip() to get a list of district ID(s).
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get officials by dristrict ID
##' \dontrun{officials <- Officials.getByDistrict(as.list(25157:25163))}
##' \dontrun{officials}
##' @export



Officials.getByDistrict <-
function (districtId) {
  
  
  # internal function
  Officials.getByDistrict.basic <- function (.districtId) {
    
    request <-  "Officials.getByDistrict?"
    inputs  <-  paste("&districtId=",.districtId,sep="")
    output  <-  pvsRequest4(request,inputs)
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(districtId, FUN= function (y) {
    
    Officials.getByDistrict.basic(.districtId=y)
    
  }
                        )
  
  
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
