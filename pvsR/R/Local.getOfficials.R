##' Get officials for a locality
##' 
##' This function is a wrapper for the Local.getOfficials() method of the PVS API Local class which returns a list of officials in a locality. The function sends a request with this method to the PVS API for all local IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Local.getOfficials(localId)
##' @param localId a character string or list of character strings with the local ID(s) (see references for details)
##' @return A data frame with a row for each official and columns with the following variables describing the official:\cr candidatelist.candidate*.candidateId,\cr candidatelist.candidate*.firstName,\cr candidatelist.candidate*.nickName,\cr candidatelist.candidate*.middleName,\cr candidatelist.candidate*.lastName,\cr candidatelist.candidate*.suffix,\cr candidatelist.candidate*.title,\cr candidatelist.candidate*.electionParties,\cr candidatelist.candidate*.electionDistrictId,\cr candidatelist.candidate*.electionStateId,\cr candidatelist.candidate*.officeParties,\cr candidatelist.candidate*.officeDistrictId,\cr candidatelist.candidate*.officeDistrictName,\cr candidatelist.candidate*.officeStateId,\cr candidatelist.candidate*.officeId,\cr candidatelist.candidate*.officeName,\cr candidatelist.candidate*.officeTypeId.
##' @references http://api.votesmart.org/docs/Local.html\cr
##' Use Local.getCounties() or Local.getCities() to get a list of local IDs. 
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of officials according to certain local IDs 
##' \dontrun{officials <- Local.getOfficials(list(3200,3203))}
##' \dontrun{officials}
##' @export


Local.getOfficials <-
function (localId) {

  
# internal function
Local.getOfficials.basic <- function (.localId) {
  
request <-  "Local.getOfficials?"
inputs  <-  paste("&localId=",.localId, sep="")
output  <-  pvsRequest4(request,inputs)
output$localId <-.localId
output

}


  # Main function
  output.list <- lapply(localId, FUN= function (b) {
    
      Local.getOfficials.basic(.localId=b)
           
        
    }
  )

output.list <- redlist(output.list)


output <- dfList(output.list)

if(class(output)=="data.frame") {
  
  output
  
}


}
