##' Get a list of candidates according to the district they represent
##' 
##' This function is a wrapper for the Candidates.getByDistrict() method of the PVS API Candidates class which grabs a list of candidates according to the district they represent for each district and election year. The function sends a request with this method to the PVS API for all district IDs and election years given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Candidates.getByDistrict(districtId, electionYear=NULL)
##' @param districtId a character string or list of character strings with the district ID(s) (see references for details)
##' @param electionYear (optional) a character string or list of character strings with the election year(s) (Default: >= current year)
##' @return A data frame with a row for each candidate and year and columns with the following variables describing the candidate:\cr candidateList.candidate*.candidateId,\cr candidateList.candidate*.firstName,\cr candidateList.candidate*.nickName,\cr candidateList.candidate*.middleName,\cr candidateList.candidate*.preferredName,\cr candidateList.candidate*.lastName,\cr candidateList.candidate*.suffix,\cr candidateList.candidate*.title,\cr candidateList.candidate*.ballotName,\cr candidateList.candidate*.electionParties,\cr candidateList.candidate*.electionStatus,\cr candidateList.candidate*.electionStage,\cr candidateList.candidate*.electionDistrictId,\cr candidateList.candidate*.electionDistrictName,\cr candidateList.candidate*.electionOffice,\cr candidateList.candidate*.electionOfficeId,\cr candidateList.candidate*.electionStateId,\cr candidateList.candidate*.electionOfficeTypeId,\cr candidateList.candidate*.electionYear,\cr candidateList.candidate*.electionSpecial,\cr candidateList.candidate*.electionDate,\cr candidateList.candidate*.officeParties,\cr candidateList.candidate*.officeStatus,\cr candidateList.candidate*.officeDistrictId,\cr candidateList.candidate*.officeDistrictName,\cr candidateList.candidate*.officeStateId,\cr candidateList.candidate*.officeId,\cr candidateList.candidate*.officeName,\cr candidateList.candidate*.officeTypeId,\cr candidateList.candidate*.runningMateId,\cr candidateList.candidate*.runningMateName.
##' @references http://api.votesmart.org/docs/Candidates.html\cr 
##' Use District.getByOfficeState() or District.getByZip() to get a list of district ID(s).
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of candidates of certain districts
##' \dontrun{district <- Candidates.getByDistrict(list(25157,25155),2012)}
##' \dontrun{district}
##' @export





Candidates.getByDistrict <-
function (districtId, electionYear=NULL) {
  
  if (length(electionYear)==0) {
    # internal function
    Candidates.getByDistrict.basic1 <- function (.districtId) {
      
      request <-  "Candidates.getByDistrict?"
      inputs  <-  paste("&districtId=",.districtId,sep="")
      output  <-  pvsRequest4(request,inputs)
      output$districtId <- .districtId
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(districtId, FUN= function (s) {
      Candidates.getByDistrict.basic1(.districtId=s)
    }
                          )
    
    
    output.list <- redlist(output.list)
    
    
    output <- dfList(output.list)
    
    
  } else {
    
    # internal function
    Candidates.getByDistrict.basic2 <- function (.districtId, .electionYear) {
      
      request <-  "Candidates.getByDistrict?"
      inputs  <-  paste("&districtId=",.districtId, "&electionYear=", .electionYear, sep="")
      output  <-  pvsRequest4(request,inputs)
      output$districtId <- .districtId
      output$electionYear.input <- .electionYear
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(districtId, FUN= function (s) {
      lapply(electionYear, FUN= function (c) {
        Candidates.getByDistrict.basic2( .districtId=s, .electionYear=c)
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
