##' Get a list of candidates according to the election they are running for
##'  
##' This function is a wrapper for the Candidates.getByElection() method of the PVS API Candidates class which grabs a list of candidates according to the election they are running for. The function sends a request with this method to the PVS API for all electionIDs and stageIDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Candidates.getByElection(electionId, stageId=NULL)
##' @param electionId a character string or list of character strings with the election ID(s) (see references for details)
##' @param stageId (optional) a character string or list of character strings with the stage ID(s) (default: all)
##' @return A data frame with a row for each candidate and columns with the following variables describing the candidate:\cr candidateList.candidate*.candidateId,\cr candidateList.candidate*.firstName,\cr candidateList.candidate*.nickName,\cr candidateList.candidate*.middleName,\cr candidateList.candidate*.preferredName,\cr candidateList.candidate*.lastName,\cr candidateList.candidate*.suffix,\cr candidateList.candidate*.title,\cr candidateList.candidate*.ballotName,\cr candidateList.candidate*.electionParties,\cr candidateList.candidate*.electionStatus,\cr candidateList.candidate*.electionStage,\cr candidateList.candidate*.electionDistrictId,\cr candidateList.candidate*.electionDistrictName,\cr candidateList.candidate*.electionOffice,\cr candidateList.candidate*.electionOfficeId,\cr candidateList.candidate*.electionStateId,\cr candidateList.candidate*.electionOfficeTypeId,\cr candidateList.candidate*.electionYear,\cr candidateList.candidate*.electionSpecial,\cr candidateList.candidate*.electionDate,\cr candidateList.candidate*.officeParties,\cr candidateList.candidate*.officeStatus,\cr candidateList.candidate*.officeDistrictId,\cr candidateList.candidate*.officeDistrictName,\cr candidateList.candidate*.officeStateId,\cr candidateList.candidate*.officeId,\cr candidateList.candidate*.officeName,\cr candidateList.candidate*.officeTypeId,\cr candidateList.candidate*.runningMateId,\cr candidateList.candidate*.runningMateName
##' @references http://api.votesmart.org/docs/Candidates.html\cr 
##' Use Election.getElectionByYearState() or Election.getElectionByZip() to get election ID(s).
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of candidates for certain election IDs
##' \dontrun{candidates <- Candidates.getByElection(list(2582,2646))}
##' \dontrun{candidates}
##' @export

Candidates.getByElection <-
function (electionId, stageId=NULL) {
  
  if (length(stageId)==0) {
    # internal function
    Candidates.getByElection.basic1 <- function (.electionId) {
      
      request <-  "Candidates.getByElection?"
      inputs  <-  paste("&electionId=",.electionId,sep="")
      output  <-  pvsRequest4(request,inputs)
      output$electionId <- .electionId
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(electionId, FUN= function (s) {
      Candidates.getByElection.basic1(.electionId=s)
    }
                          )
    
    
    output.list <- redlist(output.list)
    
    
    output <- dfList(output.list)
    
    
  } else {
    
    # internal function
    Candidates.getByElection.basic2 <- function (.electionId, .stageId) {
      
      request <-  "Candidates.getByElection?"
      inputs  <-  paste("&electionId=",.electionId, "&stageId=", .stageId, sep="")
      output  <-  pvsRequest4(request,inputs)
      output$electionId <- .electionId
      output$stageId.input <- .stageId
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(electionId, FUN= function (s) {
      lapply(stageId, FUN= function (c) {
        Candidates.getByElection.basic2( .electionId=s, .stageId=c)
      }
             )
    }
                          )
    
    
    
    output.list <- redlist(output.list)
    
    output <- dfList(output.list)
    
    
    # Avoids, that output is missleading, because stageId is already given in request-output, but also a
    # additionally generated (as stageId.input). Problem exists because some request-outputs might be empty
    # and therefore only contain one "stageId" whereas the non-empty ones contain two. (see basic function)
    output$stageId[c(as.vector(is.na(output$stageId)))] <- output$stageId.input[as.vector(is.na(output$stageId))]
    output$stageId.input <- NULL
    
    
    
  }
  
  
  
  
  
  output
  
  
  
}
