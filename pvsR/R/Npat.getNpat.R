##' Get a candidate's most recently filled out NPAT/PCT (Political Courage Test)
##' 
##' This function is a wrapper for the Npat.getNpat() method of the PVS API Npat class which returns the candidate's most recently filled out NPAT/PCT. The function sends a request with this method to the PVS API for all candidate IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Npat.getNpat(candidateId)
##' @param candidateId a character string or list of character strings with the candidate ID(s) (see references for details)
##' @return A data frame with a row for each candidate and columns with the following variables describing the candidate:\cr bio.candidate.crpId (OpenSecrets ID),\cr bio.candidate.firstName,\cr bio.candidate.nickName,\cr bio.candidate.middleName,\cr bio.candidate.lastName,\cr bio.candidate.suffix,\cr bio.candidate.birthDate,\cr bio.candidate.birthPlace,\cr bio.candidate.pronunciation,\cr bio.candidate.gender,\cr bio.candidate.family,\cr bio.candidate.photo,\cr bio.candidate.homeCity,\cr bio.candidate.homeState,\cr bio.candidate.education,\cr bio.candidate.profession,\cr bio.candidate.political,\cr bio.candidate.religion,\cr bio.candidate.congMembership,\cr bio.candidate.orgMembership,\cr bio.candidate.specialMsg,\cr bio.office.parties,\cr bio.office.title,\cr bio.office.shortTitle,\cr bio.office.name,\cr bio.office.type,\cr bio.office.status,\cr bio.office.firstElect,\cr bio.office.lastElect,\cr bio.office.nextElect,\cr bio.office.termStart,\cr bio.office.termEnd,\cr bio.office.district,\cr bio.office.districtId,\cr bio.office.stateId,\cr bio.office.committee*.committeeId,\cr bio.office.committee*.committeeName,\cr bio.election*.office,\cr bio.election*.officeId,\cr bio.election*.officeType,\cr bio.election*.parties,\cr bio.election*.district,\cr bio.election*.districtId,\cr bio.election*.status,\cr bio.election*.ballotName.
##' @references http://api.votesmart.org/docs/CandidateBio.html\cr 
##' Use Candidates.getByOfficeState(), Candidates.getByOfficeTypeState(), Candidates.getByLastname(), Candidates.getByLevenshtein(), Candidates.getByElection(), Candidates.getByDistrict() or Candidates.getByZip() to get a list of candidate IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get political courage tests of Barack Obama and John Sidney McCain III
##' \dontrun{pcts <- Npat.getNpat(list(9490,53270))}
##' \dontrun{head(pcts$survey)}
##' \dontrun{head(pcts$candidate)}
##' @export



Npat.getNpat <-
function (candidateId) {

  
# internal function 
  Npat.getNpat.basic <- function (.candidateId) {
  
request <-  "Npat.getNpat?"
inputs  <-  paste("&candidateId=",.candidateId, sep="")
output  <-  pvsRequest_PCT(request,inputs)
output  <-  lapply(output, FUN=function(i){ 
  if (length(grep(pattern="candidateId",x=names(i), value=TRUE))>0) {i}
  else {data.frame(i,"candidateId" =.candidateId)}})

}


  output.list <- lapply(candidateId, FUN= function (b) {
    
    Npat.getNpat.basic(.candidateId=b)
             
    }
  )

candidates <- lapply(output.list, function(x) x[["candidate"]])
candidates <- dfList(redlist(candidates))
surveys <- lapply(output.list, function(x) x[["survey"]])
surveys <- dfList(redlist(surveys))

return(list(candidate=candidates, survey=surveys))

}
