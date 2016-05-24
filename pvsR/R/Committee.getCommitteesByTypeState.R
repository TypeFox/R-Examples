##' Get a list of committees according to type and state
##' 
##' This function is a wrapper for the Committee.getCommitteesByTypeState() method of the PVS API Committee class which returns a list of committees for each type in each requested state. The function sends a request with this method to the PVS API for all type and state IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Committee.getCommitteesByTypeState(typeId=list("H","S","J"), stateId="NA", all=FALSE)
##' @param typeId (optional) a character string or list of character strings with the type ID(s) (default: All) (see references for details)
##' @param stateId (optional) a character string or list of character strings with the state ID(s) (default: "NA", for national) (see references for details)
##' @param all a logical indicator; if TRUE data on all possible combinations of the input variables are returned, if FALSE (default) only the exact combinations of them (see example)
##' @return A data frame with a row for each committee and columns with the following variables describing the committee:\cr committees.committee*.committeeId,\cr committees.committee*.parentId,\cr committees.committee*.stateId,\cr committees.committee*.committeeTypeId,\cr committees.committee*.name.
##' @references http://api.votesmart.org/docs/Committee.html|cr
##' See http://api.votesmart.org/docs/semi-static.html for a list of committee-type ID(s).\cr
##' Use State.getStateIDs() to get a list of state IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a data frame of committees according to all type/state combinations
##' \dontrun{committees <- Committee.getCommitteesByTypeState(typeId=list("H","S"),
##' stateId=list("NY","NJ"), all=TRUE)}
##' \dontrun{committees}
##' # get a data frame of committees according to the exact type/state combinations 
##' # (i.e., "H"/"NY", "S"/"NJ")
##' \dontrun{committees <- Committee.getCommitteesByTypeState(typeId=list("H","S"),
##' stateId=list("NY","NJ"), all=FALSE)}
##' \dontrun{committees}


##' @export



Committee.getCommitteesByTypeState <-
function (typeId=list("H","S","J"), stateId="NA", all=FALSE) {

  
# internal function
Committee.getCommitteesByTypeState.basic <- function (.typeId, .stateId) {
  
request <-  "Committee.getCommitteesByTypeState?"
inputs  <-  paste("&typeId=",.typeId,"&stateId=",.stateId,sep="")
output  <-  pvsRequest(request,inputs)
output$typeId <- .typeId
output$stateId <- .stateId
output

}  
  

if (all==TRUE) {
  

# Main function  
  output.list <- lapply(typeId, FUN= function (y) {
    lapply(stateId, FUN= function (s) {
      Committee.getCommitteesByTypeState.basic(.typeId=y, .stateId=s)
           }
        )
    }
  )

} else {
  
  # Main function  
  
  reqdf <- data.frame(t=unlist(typeId), s=unlist(stateId))
  
  output.list <- lapply(1:dim(reqdf)[1], FUN= function (l) {
    
    Committee.getCommitteesByTypeState.basic(.typeId=reqdf[l,"t"], .stateId=reqdf[l,"s"])
    
  })
  
}  

output.list <- redlist(output.list)

output <- dfList(output.list)

output


}
