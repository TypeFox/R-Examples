##' Get leadership positions by state and office
##' 
##' This function is a wrapper for the Leadership.getPositions() method of the PVS API Leadership class which returns leadership positions by state and office. The function sends a request with this method to the PVS API for all state and office IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Leadership.getPositions(stateId="NA", officeId=NULL)
##' @param stateId (optional) a character string or list of character strings with the state ID(s) (default: "NA", for national) (see references for details)
##' @param officeId (optional) a character string or list of character strings with the office ID(s) (default: All) (see references for details)
##' @return A data frame with a row for each leadership position and columns with the following variables describing the position:\cr leadership.position*.leadershipId,\cr leadership.position*.name,\cr leadership.position*.officeId,\cr leadership.position*.officeName.
##' @references http://api.votesmart.org/docs/Leadership.html\cr
##' Use State.getStateIDs() to get a list of state IDs.\cr
##' See http://api.votesmart.org/docs/semi-static.html for a list of office IDs or use Office.getOfficesByType(), Office.getOfficesByLevel(), Office.getOfficesByTypeLevel() or Office.getOfficesByBranchLevel() to get a list of office ID(s).
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get leadership positions by state ID and office ID
##' \dontrun{positions <- Leadership.getPositions(list("AL","FL"),8)}
##' \dontrun{positions}
##' @export



Leadership.getPositions <-
function (stateId="NA", officeId=NULL) {
  
  if (length(officeId)==0) {
    # internal function
    Leadership.getPositions.basic1 <- function (.stateId) {
      
      request <-  "Leadership.getPositions?"
      inputs  <-  paste("&stateId=",.stateId,sep="")
      output  <-  pvsRequest4(request,inputs)
      output$stateId <- .stateId
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(stateId, FUN= function (s) {
      Leadership.getPositions.basic1(.stateId=s)
    }
                          )
    
    
    output.list <- redlist(output.list)
    
    
    output <- dfList(output.list)
    
    
  } else {
    
    # internal function
    Leadership.getPositions.basic2 <- function (.stateId, .officeId) {
      
      request <-  "Leadership.getPositions?"
      inputs  <-  paste("&stateId=",.stateId, "&officeId=", .officeId, sep="")
      output  <-  pvsRequest4(request,inputs)
      output$stateId <- .stateId
      output$officeId.input <- .officeId
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(stateId, FUN= function (s) {
      lapply(officeId, FUN= function (c) {
        Leadership.getPositions.basic2( .stateId=s, .officeId=c)
      }
             )
    }
                          )
    
    
    
    output.list <- redlist(output.list)
    
    output <- dfList(output.list)
    
    
    # Avoids, that output is missleading, because officeId is already given in request-output, but also a
    # additionally generated (as officeId.input). Problem exists because some request-outputs might be empty
    # and therefore only contain one "officeId" whereas the non-empty ones contain two. (see basic function)
    output$officeId[c(as.vector(is.na(output$officeId)))] <- output$officeId.input[as.vector(is.na(output$officeId))]
    output$officeId.input <- NULL
    
    
    
  }
  
  
  
  
  
  output
  
  
  
}
