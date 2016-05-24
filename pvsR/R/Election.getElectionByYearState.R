##' Get district basic election data according to year and state ID
##' 
##' This function is a wrapper for the Election.getElectionByYearState() method of the PVS API Election class which grabs district basic election data. The function sends a request with this method to the PVS API for all state IDs and years given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Election.getElectionByYearState(stateId="NA", year)
##' @param stateId (optional) a character string or list of character strings with the state ID(s) (default: "NA", for national) (see references for details)
##' @param year a character string or list of character strings with the years (defaults to current year)
##' @return A data frame with a row for each election and columns with the following variables describing the election:\cr elections.election*.electionId,\cr elections.election*.name,\cr elections.election*.stateId,\cr elections.election*.officeTypeId,\cr elections.election*.special,\cr elections.election*.electionYear.\cr For each stage the following variables are jointly (as one string) in a column:\cr elections.election*.stage*.stageId,\cr elections.election*.stage*.name,\cr elections.election*.stage*.stateId,\cr elections.election*.stage*.electionDate,\cr elections.election*.stage*.filingDeadline,\cr elections.election*.stage*.npatMailed.
##' @references http://api.votesmart.org/docs/Election.html\cr
##' Use State.getStateIDs() to get a list of state IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get information about an election of a certain year and state 
##' \dontrun{election <- Election.getElectionByYearState(list("NY","FL"),2012)}
##' \dontrun{election}
##' @export




Election.getElectionByYearState <-
function (stateId="NA", year=NULL) {
  
  if (is.null(year)) year <- substr(Sys.Date(), 1,4)
  
  # internal function
  Election.getElectionByYearState.basic <- function (.stateId, .year) {
    
    request <-  "Election.getElectionByYearState?"
    inputs  <-  paste("&stateId=",.stateId,"&year=",.year,sep="")
    output  <-  pvsRequest9.1(request,inputs)
    output$stateId <- .stateId
    output$year <- .year
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(stateId, FUN= function (y) {
    lapply(year, FUN= function (s) {
      Election.getElectionByYearState.basic(.stateId=y, .year=s)
    }
           )
  }
                        )
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
