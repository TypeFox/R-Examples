##' Get district basic election data according to the ZIP code
##' 
##' This function is a wrapper for the Election.getElectionByZip() method of the PVS API Election class which grabs district basic election data according to the ZIP code. If another year than the current year is chosen, all election data from that year up to the current year is returned. The function sends a request with this method to the PVS API for all ZIP codes given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Election.getElectionByZip(zip5, zip4=NULL, year=NULL)
##' @param zip5 a character string or list of character strings with the five-digit ZIP code
##' @param zip4 (optional) a character string or list of character strings with the expanded ZIP+4 code (default: NULL)
##' @param year a character string or list of character strings with the year (defaults to current year)
##' @return A data frame with a row for each election and columns with the following variables describing the election:\cr elections.election*.electionId,\cr elections.election*.name,\cr elections.election*.stateId,\cr elections.election*.officeTypeId,\cr elections.election*.special,\cr elections.election*.electionYear.
##' @references http://api.votesmart.org/docs/Election.html
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get election data by ZIP code 
##' \dontrun{election <- Election.getElectionByZip(zip5=list(10001,10002), year="2012")}
##' \dontrun{election}
##' @export




Election.getElectionByZip <-
function (zip5, zip4=NULL, year=NULL) {
  
  if (is.null(year)) year <- substr(Sys.Date(), 1,4)
  
  
  if (length(zip4)==0) {
    # internal function
    Election.getElectionByZip.basic1 <- function (.zip5, .year) {
      
      request <-  "Election.getElectionByZip?"
      inputs  <-  paste("&zip5=",.zip5, "&year=", .year, sep="")
      output  <-  pvsRequest6(request,inputs)
      output$zip5 <- .zip5
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(zip5, FUN= function (s) {
      lapply(year, function(y) {
      Election.getElectionByZip.basic1(.zip5=s, .year=y)
      }
      )
    }
                          )
    
    
    output.list <- redlist(output.list)
    
    
    output <- dfList(output.list)
    
    
  } else {
    
    # internal function
    Election.getElectionByZip.basic2 <- function (.zip5, .zip4, .year) {
      
      request <-  "Election.getElectionByZip?"
      inputs  <-  paste("&zip5=",.zip5, "&zip4=", .zip4, "&year=", .year, sep="")
      output  <-  pvsRequest5(request,inputs)
      output$zip5 <- .zip5
      output$zip4.input <- .zip4
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(zip5, FUN= function (s) {
      lapply(zip4, FUN= function (c) {
        lapply(year, FUN= function(y) { 
        Election.getElectionByZip.basic2( .zip5=s, .zip4=c, .year=y)
        }
        )
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
