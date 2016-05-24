##' Get a list of states and their IDs.
##' 
##' This function is a wrapper for the State.getStateIDs() method of the PVS API State class which returns a simple state ID and name list for mapping IDs to state names. 
##' @usage State.getStateIDs()
##' @return A data frame with a row for each state:\cr
##' statelist.list.state*.stateId,\cr
##' statelist.list.state*.name
##' @references http://api.votesmart.org/docs/State.html
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of states with their IDs
##' \dontrun{stateIDs <- State.getStateIDs()}
##' \dontrun{stateIDs}
##' @export



State.getStateIDs <-
function () {
  
  
  
  
  states <- xmlTreeParse(paste("http://api.votesmart.org/State.getStateIDs?key=",get('pvs.key',envir=.GlobalEnv),sep=""), useInternalNodes = TRUE)
  states <- xmlRoot(states)
  states <- states[["list"]]
  
  
  statelist <- xmlSApply(states, function(x) xmlSApply(x, xmlValue))
  
  
  states.df <- as.data.frame(t(statelist),row.names=FALSE)
  
  
  states.df
  
  
}
