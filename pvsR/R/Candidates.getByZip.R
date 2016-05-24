##' Get a list of candidates according to ZIP code
##' 
##' This function is a wrapper for the Candidates.getByZip() method of the PVS API Candidates class which grabs a list of candidates according to the ZIP code. The function sends a request with this method to the PVS API for all ZIP codes and election years given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Candidates.getByZip(zip5, electionYear=NULL)
##' @param zip5 a character string or list of character strings with the 5-digit ZIP code(s)
##' @param electionYear (optional) a character string or list of character strings with the election year(s) (default: >= current year)
##' @return A data frame with a row for each candidate and columns with the following variables describing the candidate:\cr candidateList.candidate*.candidateId,\cr candidateList.candidate*.firstName,\cr candidateList.candidate*.nickName,\cr candidateList.candidate*.middleName,\cr candidateList.candidate*.preferredName,\cr candidateList.candidate*.lastName,\cr candidateList.candidate*.suffix,\cr candidateList.candidate*.title,\cr candidateList.candidate*.ballotName,\cr candidateList.candidate*.electionParties,\cr candidateList.candidate*.electionStatus,\cr candidateList.candidate*.electionStage,\cr candidateList.candidate*.electionDistrictId,\cr candidateList.candidate*.electionDistrictName,\cr candidateList.candidate*.electionOffice,\cr candidateList.candidate*.electionOfficeId,\cr candidateList.candidate*.electionStateId,\cr candidateList.candidate*.electionOfficeTypeId,\cr candidateList.candidate*.electionYear,\cr candidateList.candidate*.electionSpecial,\cr candidateList.candidate*.electionDate,\cr candidateList.candidate*.officeParties,\cr candidateList.candidate*.officeStatus,\cr candidateList.candidate*.officeDistrictId,\cr candidateList.candidate*.officeDistrictName,\cr candidateList.candidate*.officeStateId,\cr candidateList.candidate*.officeId,\cr candidateList.candidate*.officeName,\cr candidateList.candidate*.officeTypeId,\cr candidateList.candidate*.runningMateId,\cr candidateList.candidate*.runningMateName.
##' @references http://api.votesmart.org/docs/Candidates.html
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of candidates according to the ZIP code.
##' \dontrun{candidates <- Candidates.getByZip(list(10001,10002),2012)}
##' \dontrun{candidates}

##' @export



Candidates.getByZip <-
function (zip5, electionYear=NULL) {
  
  if (length(electionYear)==0) {
    # internal function
    Candidates.getByZip.basic1 <- function (.zip5) {
      
      request <-  "Candidates.getByZip?"
      inputs  <-  paste("&zip5=",.zip5,sep="")
      output  <-  pvsRequest6(request,inputs)
      output$zip5 <- .zip5
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(zip5, FUN= function (s) {
      Candidates.getByZip.basic1(.zip5=s)
    }
                          )
    
    
    output.list <- redlist(output.list)
    
    
    output <- dfList(output.list)
    
    
  } else {
    
    # internal function
    Candidates.getByZip.basic2 <- function (.zip5, .electionYear) {
      
      request <-  "Candidates.getByZip?"
      inputs  <-  paste("&zip5=",.zip5, "&electionYear=", .electionYear, sep="")
      output  <-  pvsRequest6(request,inputs)
      output$zip5 <- .zip5
      output$electionYear.input <- .electionYear
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(zip5, FUN= function (s) {
      lapply(electionYear, FUN= function (y) {
        Candidates.getByZip.basic2( .zip5=s, .electionYear=y)
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
