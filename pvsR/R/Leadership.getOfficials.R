##' Get officials that hold the leadership role in certain states
##' 
##' This function is a wrapper for the Leadership.getOfficials() method of the PVS API Leadership class which grabs a list of officials that hold the leadership role in certain states. The function sends a request with this method to the PVS API for all state and leadership IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Leadership.getOfficials(stateId="NA", leadershipId)
##' @param stateId (optional) a character string or list of character strings with the state ID(s) (default: "NA", for national) (see references for details)
##' @param leadershipId a character string or list of character strings with the leadership ID(s) (see references for details)
##' @return A data frame with a row for each leadership position and columns with the following variables describing the official:\cr leaders.leader*.candidateId,\cr leaders.leader*.firstName,\cr leaders.leader*.middleName,\cr leaders.leader*.lastName,\cr leaders.leader*.suffix,\cr leaders.leader*.position,\cr leaders.leader*.officeId,\cr leaders.leader*.title.
##' @references http://api.votesmart.org/docs/Leadership.html\cr
##' Use State.getStateIDs() to get a list of state IDs.\cr
##' Use Leadership.getPositions() to get a list of leadership IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get leaders by state ID and leadership ID 
##' \dontrun{officials <- Leadership.getOfficials(list("NY","FL"),list(138,140))}
##' \dontrun{officials}
##' @export


Leadership.getOfficials <-
function (stateId="NA", leadershipId) {
  
  
  # internal function
  Leadership.getOfficials.basic <- function (.stateId, .leadershipId) {
    
    request <-  "Leadership.getOfficials?"
    inputs  <-  paste("&stateId=",.stateId,"&leadershipId=",.leadershipId,sep="")
    output  <-  pvsRequest4(request,inputs)
    output$stateId <- .stateId
    output$leadershipId <- .leadershipId
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(stateId, FUN= function (y) {
    lapply(leadershipId, FUN= function (s) {
      Leadership.getOfficials.basic(.stateId=y, .leadershipId=s)
    }
           )
    
  }
                        
                        )
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
