##' Get district basic election data according to election ID and stage ID
##' 
##' This function is a wrapper for the Election.getStageCandidates() method of the PVS API Election class which grabs district basic election data according to the election and stage ID. The function sends a request with this method to the PVS API for all election and stage IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Election.getStageCandidates(stageId, electionId)
##' @param electionId a character string or list of character strings with the election ID(s) (see references for details)
##' @param stageId a character string or list of character strings with the stage ID(s)
##' @return A data frame with a row for each combination of stage ID and election and columns with the following variables describing the election:\cr stagecandidates.election*.electionId,\cr stagecandidates.election*.name,\cr stagecandidates.election*.stage,\cr stagecandidates.election*.stateId,\cr stagecandidates.candidate*.candidateId,\cr stagecandidates.candidate*.officeId,\cr stagecandidates.candidate*.district,\cr stagecandidates.candidate*.firstName,\cr stagecandidates.candidate*.middleName,\cr stagecandidates.candidate*.lastName,\cr stagecandidates.candidate*.suffix,\cr stagecandidates.candidate*.party,\cr stagecandidates.candidate*.status,\cr stagecandidates.candidate*.voteCount,\cr stagecandidates.candidate*.votePercent.
##' @references http://api.votesmart.org/docs/Election.html\cr
##' Use Election.getElectionByYearState() or Election.getElectionByZip() to get election ID(s).
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get election data for certain election IDs
##' \dontrun{election <- Election.getStageCandidates(stageId=list("P"),electionId=list(2582,2607))}
##' \dontrun{election}
##' @export



Election.getStageCandidates <-
function (stageId, electionId) {
  
  
  # internal function
  Election.getStageCandidates.basic <- function (.stageId, .electionId) {
    
    request <-  "Election.getStageCandidates?"
    inputs  <-  paste("&stageId=",.stageId,"&electionId=",.electionId,sep="")
    output  <-  pvsRequest6(request,inputs)
    output$stageId <- .stageId
    output$electionId <- .electionId
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(stageId, FUN= function (y) {
    lapply(electionId, FUN= function (s) {
      Election.getStageCandidates.basic(.stageId=y, .electionId=s)
    }
    )

                        }

                        )
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
