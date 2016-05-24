##' Get a list of candidates according to their lastname
##'  
##' This function is a wrapper for the Candidates.getByLastname() method of the PVS API Candidates class which grabs a list of candidates according to a lastname match. The function sends a request with this method to the PVS API for all last names and election years given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Candidates.getByLastname(lastName, electionYear=NULL)
##' @param lastName a character string or list of character strings with the last name(s) (see references for details)
##' @param electionYear (optional) a character string or list of character strings with the election year(s) (default: >= current year)
##' @return A data frame with a row for each candidate and year and columns with the following variables describing the candidate:\cr candidateList.candidate*.candidateId,\cr candidateList.candidate*.firstName,\cr candidateList.candidate*.nickName,\cr candidateList.candidate*.middleName,\cr candidateList.candidate*.preferredName,\cr candidateList.candidate*.lastName,\cr candidateList.candidate*.suffix,\cr candidateList.candidate*.title,\cr candidateList.candidate*.ballotName,\cr candidateList.candidate*.electionParties,\cr candidateList.candidate*.electionStatus,\cr candidateList.candidate*.electionStage,\cr candidateList.candidate*.electionDistrictId,\cr candidateList.candidate*.electionDistrictName,\cr candidateList.candidate*.electionOffice,\cr candidateList.candidate*.electionOfficeId,\cr candidateList.candidate*.electionStateId,\cr candidateList.candidate*.electionOfficeTypeId,\cr candidateList.candidate*.electionYear,\cr candidateList.candidate*.electionSpecial,\cr candidateList.candidate*.electionDate,\cr candidateList.candidate*.officeParties,\cr candidateList.candidate*.officeStatus,\cr candidateList.candidate*.officeDistrictId,\cr candidateList.candidate*.officeDistrictName,\cr candidateList.candidate*.officeStateId,\cr candidateList.candidate*.officeId,\cr candidateList.candidate*.officeName,\cr candidateList.candidate*.officeTypeId,\cr candidateList.candidate*.runningMateId,\cr candidateList.candidate*.runningMateName.
##' @references http://api.votesmart.org/docs/Candidates.html\cr 
##' Use CandidateBio.getBio(), Candidates.getByOfficeState(), Candidates.getByOfficeTypeState(), Candidates.getByElection(), Candidates.getByDistrict(), Candidates.getByZip(), Committee.getCommitteeMembers(), Election.getStageCandidates(), Leadership.getOfficials(), Local.getOfficials(), Officials.getStatewide(), Officials.getByOfficeState(), Officials.getByOfficeTypeState(), Officials.getByDistrict() or Officials.getByZip() to get last name(s).
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get candidates by last names
##' \dontrun{names <- Candidates.getByLastname(list("Obama","Romney"),2012)}
##' \dontrun{names}
##' @export



Candidates.getByLastname <-
function (lastName, electionYear=NULL) {
  
  if (length(electionYear)==0) {
    # internal function
    Candidates.getByLastname.basic1 <- function (.lastName) {
      
      request <-  "Candidates.getByLastname?"
      inputs  <-  paste("&lastName=",.lastName,sep="")
      output  <-  pvsRequest4(request,inputs)
      output$lastName <- .lastName
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(lastName, FUN= function (s) {
      Candidates.getByLastname.basic1(.lastName=s)
    }
                          )
    
    
    output.list <- redlist(output.list)
    
    
    output <- dfList(output.list)
    
    
  } else {
    
    # internal function
    Candidates.getByLastname.basic2 <- function (.lastName, .electionYear) {
      
      request <-  "Candidates.getByLastname?"
      inputs  <-  paste("&lastName=",.lastName, "&electionYear=", .electionYear, sep="")
      output  <-  pvsRequest4(request,inputs)
      output$lastName <- .lastName
      output$electionYear.input <- .electionYear
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(lastName, FUN= function (s) {
      lapply(electionYear, FUN= function (c) {
        Candidates.getByLastname.basic2( .lastName=s, .electionYear=c)
      }
             )
    }
                          )
    
    
    
    output.list <- redlist(output.list)
    
    output <- dfList(output.list)
    
    
    # Avoids, that output is missleading, because electionYear is already given in request-output, but also a
    # additionally generated (as electionYear.input). Problem exists because some request-outputs might be empty
    # and therefore only contain one "electionYear" whereas the non-empty ones contain two. (see basic function)
    output$electionYear[c(as.vector(is.na(output$electionYear)))] <- output$electionYear.input[as.vector(is.na(output$electionYear))]
    output$electionYear.input <- NULL
    
    
    
  }
  
  
  
  
  
  output
  
  
  
}
