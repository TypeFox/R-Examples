##' Get information about a state
##' 
##' This function is a wrapper for the State.getState() method of the PVS API State class which grabs various data on a state. The function sends a request with this method to the PVS API for all state IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage State.getState(stateId)
##' @param stateId a character string or list of character strings with the state ID(s) (see references for details)
##' @return A data frame with a row for each state and columns with the following variables describing the state:\cr state.details.stateId,\cr state.details.stateType,\cr state.details.name,\cr state.details.nickName,\cr state.details.capital,\cr state.details.area,\cr state.details.population,\cr state.details.statehood,\cr state.details.motto,\cr state.details.flower,\cr state.details.tree,\cr state.details.bird,\cr state.details.highPoint,\cr state.details.lowPoint,\cr state.details.bicameral,\cr state.details.upperLegis,\cr state.details.lowerLegis,\cr state.details.ltGov,\cr state.details.senators,\cr state.details.reps,\cr state.details.termLimit,\cr state.details.termLength,\cr state.details.billUrl,\cr state.details.voteUrl,\cr state.details.voterReg,\cr state.details.primaryDate,\cr state.details.generalDate,\cr state.details.absenteeWho,\cr state.details.absenteeHow,\cr state.details.absenteeWhen,\cr state.details.largestCity,\cr state.details.rollUpper,\cr state.details.rollLower,\cr state.details.usCircuit.
##' @references http://api.votesmart.org/docs/State.html\cr
##' Use State.getStateIDs() to get a list of state IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get information about certain states
##' \dontrun{stateinfo <- State.getState(list("FL","NY"))}
##' \dontrun{stateinfo}

##' @export




State.getState <-
function (stateId) {
  
  
  # internal function
  State.getState.basic <- function (.stateId) {
    
    request <-  "State.getState?"
    inputs  <-  paste("&stateId=",.stateId,sep="")
    output  <-  pvsRequest4(request,inputs)
    output$stateId <- .stateId
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(stateId, FUN= function (s) {
    State.getState.basic(.stateId=s)
  }
                        )
  
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
