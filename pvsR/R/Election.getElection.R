##' Get district basic election data
##' 
##' This function is a wrapper for the Election.getElection() method of the PVS API Election class which grabs district basic election data according to the election ID. The function sends a request with this method to the PVS API for all election IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Election.getElection(electionId)
##' @param electionId a character string or list of character strings with the election ID(s) (see references for details)
##' @return A data frame with a row for each election and columns with the following variables describing the election:\cr elections.election*.electionId,\cr elections.election*.name,\cr elections.election*.stateId,\cr elections.election*.officeTypeId,\cr elections.election*.special,\cr elections.election*.electionYear.\cr For each stage the following variables are jointly (as one string) in a column:\cr elections.election*.stage*.stageId,\cr elections.election*.stage*.name,\cr elections.election*.stage*.stateId,\cr elections.election*.stage*.electionDate,\cr elections.election*.stage*.filingDeadline,\cr elections.election*.stage*.npatMailed.
##' @references http://api.votesmart.org/docs/Election.html\cr
##' Use Election.getElectionByYearState() or Election.getElectionByZip() to get election ID(s).
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get data of a certain election
##' \dontrun{election <- Election.getElection(as.list(2595:2607))}
##' \dontrun{election}
##' @export





Election.getElection <-
function (electionId) {
  
  
  # internal function
  Election.getElection.basic <- function (.electionId) {
    
    request <-  "Election.getElection?"
    inputs  <-  paste("&electionId=",.electionId,sep="")
    output  <-  pvsRequest10(request,inputs)
    
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(electionId, FUN= function (s) {
    Election.getElection.basic(.electionId=s)
  }
                        )
  
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
