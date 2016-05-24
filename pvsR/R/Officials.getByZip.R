##' Get a list of officials according to the ZIP code
##' 
##' This function is a wrapper for the Officials.getByZip() method of the PVS API Officials class which grabs a list of officials according to the ZIP code of the area they represent. The function sends a request with this method to the PVS API for all ZIP Codes given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Officials.getByZip(zip5, zip4=NULL)
##' @param zip5 a character string or list of character strings with the five-digit ZIP code
##' @param zip4 (optional) a character string or list of character strings with the expanded ZIP+4 code (default: all)
##' @return A data frame with a row for each official and columns with the following variables describing the official:\cr candidateList.zipMessage,\cr candidateList.candidate*.candidateId,\cr candidateList.candidate*.firstName,\cr candidateList.candidate*.nickName,\cr candidateList.candidate*.middleName,\cr candidateList.candidate*.lastName,\cr candidateList.candidate*.suffix,\cr candidateList.candidate*.title,\cr candidateList.candidate*.electionParties,\cr candidateList.candidate*.electionstatus,\cr candidateList.candidate*.officeParties,\cr candidatelist.candidate*.officeStatus,\cr candidateList.candidate*.officeDistrictId,\cr candidateList.candidate*.officeDistrictName,\cr candidateList.candidate*.officeTypeId,\cr candidateList.candidate*.officeId,\cr candidateList.candidate*.officeName,\cr candidateList.candidate*.officeStateId.
##' @references http://api.votesmart.org/docs/Officials.html
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of officials by ZIP code
##' \dontrun{officials <- Officials.getByZip(list(10001,10002))}
##' \dontrun{head(officials)}

##' @export







Officials.getByZip <-
function (zip5, zip4=NULL) {
  
  if (length(zip4)==0) {
    # internal function
    Officials.getByZip.basic1 <- function (.zip5) {
      
      request <-  "Officials.getByZip?"
      inputs  <-  paste("&zip5=",.zip5,sep="")
      output  <-  pvsRequest5(request,inputs)
      output$zip5 <- .zip5
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(zip5, FUN= function (s) {
      Officials.getByZip.basic1(.zip5=s)
    }
                          )
    
    
    output.list <- redlist(output.list)
    
    
    output <- dfList(output.list)
    
    
  } else {
    
    # internal function
    Officials.getByZip.basic2 <- function (.zip5, .zip4) {
      
      request <-  "Officials.getByZip?"
      inputs  <-  paste("&zip5=",.zip5, "&zip4=", .zip4, sep="")
      output  <-  pvsRequest4(request,inputs)
      output$zip5 <- .zip5
      output$zip4.input <- .zip4
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(zip5, FUN= function (s) {
      lapply(zip4, FUN= function (c) {
        Officials.getByZip.basic2( .zip5=s, .zip4=c)
      }
             )
    }
                          )
    
    
    
    output.list <- redlist(output.list)
    
    output <- dfList(output.list)
    
    
    # Avoids, that output is missleading, because zip4 is already given in request-output, but also a
    # additionally generated (as zip4.input). Problem exists because some request-outputs might be empty
    # and therefore only contain one "zip4" whereas the non-empty ones contain two. (see basic function)
    output$zip4[c(as.vector(is.na(output$zip4)))] <- output$zip4.input[as.vector(is.na(output$zip4))]
    output$zip4.input <- NULL
    
    
    
  }
  
  
  
  
  
  output
  
  
  
}
