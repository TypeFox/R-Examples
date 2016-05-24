##' Get a list of officials according to a last name match
##' 
##' This function is a wrapper for the Officials.getByLastname() method of the PVS API Officials class which grabs a list of officials according to a last name match. The function sends a request with this method to the PVS API for all last names given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Officials.getByLastname(lastName)
##' @param lastName a character string or list of character strings with the last name(s) (see references for details)
##' @return A data frame with a row for each official and columns with the following variables describing the official:\cr candidateList.candidate*.candidateId,\cr candidateList.candidate*.firstName,\cr candidateList.candidate*.nickName,\cr candidateList.candidate*.middleName,\cr candidateList.candidate*.lastName,\cr candidateList.candidate*.suffix,\cr candidateList.candidate*.title,\cr candidateList.candidate*.electionParties,\cr candidateList.candidate*.officeParties,\cr candidatelist.candidate*.officeStatus,\cr candidateList.candidate*.officeDistrictId,\cr candidateList.candidate*.officeDistrictName,\cr candidateList.candidate*.officeTypeId,\cr candidateList.candidate*.officeId,\cr candidateList.candidate*.officeName,\cr candidateList.candidate*.officeStateId.
##' @references http://api.votesmart.org/docs/Officials.html\cr
##' Use CandidateBio.getBio(), Candidates.getByOfficeState(), Candidates.getByOfficeTypeState(), Candidates.getByElection(), Candidates.getByDistrict(), Candidates.getByZip()\cr, Committee.getCommitteeMembers(), Election.getStageCandidates(), Leadership.getOfficials(), Local.getOfficials(), Officials.getStatewide(), Officials.getByOfficeState(), Officials.getByOfficeTypeState(), Officials.getByDistrict() or Officials.getByZip() to get last name(s).
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of officials with the same last name
##' \dontrun{miller <- Officials.getByLastname(list("Miller","Fine"))}
##' \dontrun{miller}
##' @export



Officials.getByLastname <-
function (lastName) {
  
  
  # internal function
  Officials.getByLastname.basic <- function (.lastName) {
    
    request <-  "Officials.getByLastname?"
    inputs  <-  paste("&lastName=",.lastName,sep="")
    output  <-  pvsRequest4(request,inputs)
    output$lastName <- .lastName
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(lastName, FUN= function (y) {
    
    Officials.getByLastname.basic(.lastName=y)
    
  }
                        )
  
  
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
