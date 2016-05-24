##' Get counties in a state
##' 
##' This function is a wrapper for the Local.getCounties() method of the PVS API Local class which returns a list of counties in a state. The function sends a request with this method to the PVS API for all state IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Local.getCounties(stateId)
##' @param stateId a character string or list of character strings with the state ID(s) (see references for details)
##' @return A data frame with a row for each county and columns with the following variables describing the county:\cr counties.county*.localId,\cr counties.county*.name,\cr counties.county*.url.
##' @references http://api.votesmart.org/docs/Local.html\cr
##' Use State.getStateIDs() to get a list of state IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of counties of a certain state
##' \dontrun{counties <- Local.getCounties(list("NY","FL"))}
##' \dontrun{head(counties)}
##' @export



Local.getCounties <-
function (stateId) {

  
# internal function
Local.getCounties.basic <- function (.stateId) {
  
request <-  "Local.getCounties?"
inputs  <-  paste("&stateId=",.stateId, sep="")
output  <-  pvsRequest4(request,inputs)
output$stateId <-.stateId
output

}


  # Main function
  output.list <- lapply(stateId, FUN= function (b) {
    
      Local.getCounties.basic(.stateId=b)
           
        
    }
  )

output.list <- redlist(output.list)


output <- dfList(output.list)

if(class(output)=="data.frame") {
  
  output
  
}




}
