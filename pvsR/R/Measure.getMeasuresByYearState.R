##' Get a list of state ballot measures in a given year
##' 
##' This function is a wrapper for the Measure.getMeasuresByYearState() method of the PVS API Measure class which grabs a list of state ballot measures in a given year. The function sends a request with this method to the PVS API for all years and state IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Measure.getMeasuresByYearState(year, stateId)
##' @param year a character string or list of character strings with the year(s)
##' @param stateId a character string or list of character strings with the state ID(s) (see references for details)
##' @return A data frame with a row for each ballot measure and columns with the following variables describing the ballot measure:\cr measures.measure*.measureId,\cr measures.measure*.measureCode,\cr measures.measure*.title,\cr measures.measure*.outcome.
##' @references http://api.v.otesmart.org/docs/Measure.html\cr
##' Use State.getStateIDs() to get a list of state IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of ballot measures for a certain state and year
##' \dontrun{measures <- Measure.getMeasuresByYearState(list(2010,2012),"FL")}
##' \dontrun{measures}

##' @export






Measure.getMeasuresByYearState <-
function (year, stateId) {
  

# internal function
Measure.getMeasuresByYearState.basic <- function (.year, .stateId) {
  
request <-  "Measure.getMeasuresByYearState?"
inputs  <-  paste("&year=",.year,"&stateId=",.stateId,sep="")
output  <-  pvsRequest(request,inputs)
output$year <-.year
output$stateId <- .stateId
output

}


  #Main function
  output.list <- lapply(year, FUN= function (y) {
    lapply(stateId, FUN= function (s) {
      Measure.getMeasuresByYearState.basic(.year=y, .stateId=s)
           }
        )
    }
  )


output.list <- do.call("c",output.list)


# which list entry has the most columns, how many are these?
coln <- which.is.max(sapply(output.list, ncol));
max.cols <- max(sapply(output.list, ncol));

# give all list entries (dfs in list) the same number of columns and the same names
output.list2 <- lapply(output.list, function(x){
if (ncol(x) < max.cols) x <- data.frame(cbind(matrix(NA, ncol=max.cols-ncol(x), nrow = 1, ),x),row.names=NULL)
names(x) <- names(output.list[[coln]])
x
})

output <- do.call("rbind",output.list2)
output


}
