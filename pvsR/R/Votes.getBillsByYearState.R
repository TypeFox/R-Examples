##' Get a list of bills according to year and state
##' 
##' This function is a wrapper for the Votes.getBillsByYearState() method of the PVS API Votes class which returns a list of bills that fit the year and state input. The function sends a request with this method to the PVS API for all years and state IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Votes.getBillsByYearState(year, stateId, all=FALSE)
##' @param year a character string or list of character strings with the year(s)
##' @param stateId a character string or list of character strings with the state ID(s) (see references for details)
##' @param all a logical indicator; if TRUE data on all possible combinations of the stateId and year are returned, if FALSE (default) only the exact combinations (see example)
##' @return A data frame with a row for each bill and columns with the following variables describing the bill:\cr bills.bill*.billId,\cr bills.bill*.billNumber,\cr bills.bill*.title,\cr bills.bill*.type.
##' @references http://api.votesmart.org/docs/Votes.html\cr
##' Use State.getStateIDs() to get a list of state IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a data frame of bills according to all year and state combinations
##' \dontrun{bills <- Votes.getBillsByYearState(year=list(2011,2012),
##' stateId=list("NY","NJ"), all=TRUE)}
##' \dontrun{head(bills)}
##' # get a data frame of bills according to the exact year and state combinations 
##' # (i.e., 2011/"NY", 2012/"NJ")
##' \dontrun{bills <- Votes.getBillsByYearState(year=list(2011,2012),
##' stateId=list("NY","NJ"), all=FALSE)}
##' \dontrun{head(bills)}
##' @export



Votes.getBillsByYearState <-
function (year, stateId, all=FALSE) {
  
 
    
# internal function
Votes.getBillsByYearState.basic <- function (.year, .stateId) {
  
request <-  "Votes.getBillsByYearState?"
inputs  <-  paste("&year=",.year,"&stateId=",.stateId,sep="")
output  <-  pvsRequest4(request,inputs)
output$year <- .year
output$stateId <- .stateId
output

}
  
if (all==TRUE) {

# Main function  

  output.list <- lapply(year, FUN= function (y) {
    lapply(stateId, FUN= function (s) {
      Votes.getBillsByYearState.basic(.year=y, .stateId=s)
           }
        )
    }
  )


} else {

  # Main function  
  
  reqdf <- data.frame(y=unlist(year), s=unlist(stateId))
  
  output.list <- lapply(1:dim(reqdf)[1], FUN= function (l) {
    
      Votes.getBillsByYearState.basic(.year=reqdf[l,"y"], .stateId=reqdf[l,"s"])

  })

}  
  
output.list <- redlist(output.list)

output <- dfList(output.list)

output


}
