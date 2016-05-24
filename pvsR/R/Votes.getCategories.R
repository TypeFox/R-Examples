##' Get a list of categories that contain released bills according to year and state
##' 
##' This function is a wrapper for the Votes.getCategories() method of the PVS API Votes class which dumps categories that contain released bills according to year and state. The function sends a request with this method to the PVS API for all years and state IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Votes.getCategories(year, stateId="NA")
##' @param year a character string or list of character strings with the year(s)
##' @param stateId (optional) a character string or list of character strings with the state ID(s) (default: "NA", for national) (see references for details)
##' @return A data frame with a row for each category and columns with the following variables describing the category:\cr categories.category*.categoryId,\cr categories.category*.name.
##' @references http://api.votesmart.org/docs/Votes.html\cr
##' Use State.getStateIDs() to get a list of state IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of categories
##' \dontrun{categories <- Votes.getCategories(2012)}
##' \dontrun{categories}

##' @export




Votes.getCategories <-
function (year, stateId="NA") {
  

# internal function
Votes.getCategories.basic <- function (.year, .stateId) {
  
request <-  "Votes.getCategories?"
inputs  <-  paste("&year=",.year,"&stateId=",.stateId,sep="")
output  <-  pvsRequest(request,inputs)
output$year <-.year
output$stateId <- .stateId
output

}


  #Main function
  output.list <- lapply(year, FUN= function (y) {
    lapply(stateId, FUN= function (s) {
      Votes.getCategories.basic(.year=y, .stateId=s)
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
