##' Get a list of bills according to category, year and state
##' 
##' This function is a wrapper for the Votes.getBillsByCategoryYearState() method of the PVS API Votes class which grabs a list of bills according to category, year and state. The function sends a request with this method to the PVS API for all years, state and category IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Votes.getBillsByCategoryYearState(year, stateId, categoryId)
##' @param year a character string or list of character strings with the year ID(s)
##' @param stateId a character string or list of character strings with the state ID(s) (see references for details)
##' @param categoryId a character string or list of character strings with the category ID(s) (see references for details)
##' @return A data frame with a row for each bill and columns with the following variables describing the bill:\cr bills.bill*.billId,\cr bills.bill*.billNumber,\cr bills.bill*.title,\cr bills.bill*.type.
##' @references http://api.votesmart.org/docs/Votes.html\cr
##' Use State.getStateIDs() to get a list of state IDs.\cr
##' Use Votes.getCategories() or Rating.getCandidateRating() to get a list of category IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of bills in a certain category, year and state
##' \dontrun{bills <- Votes.getBillsByCategoryYearState(as.list(2010:2012),"NY",10)}
##' \dontrun{bills}
##' @export



Votes.getBillsByCategoryYearState <-
function (year, stateId, categoryId) {
  
    
# internal function
Votes.getBillsByCategoryYearState.basic <- function (.year, .stateId, .categoryId) {
  
request <-  "Votes.getBillsByCategoryYearState?"
inputs  <-  paste("&year=",.year,"&stateId=",.stateId, "&categoryId=", .categoryId, sep="")
output  <-  pvsRequest(request,inputs)
output$year <- .year
output$stateId <- .stateId
output$categoryId <- .categoryId
output

}
  

# Main function  

  output.list <- lapply(year, FUN= function (y) {
    lapply(stateId, FUN= function (s) {
      lapply(categoryId, FUN= function (c) {
       Votes.getBillsByCategoryYearState.basic(.year=y, .stateId=s, .categoryId=c)
              }
             )
           }
        )
    }
  )
  
  
# reduce lists in list (3fold) to one list
output.list <- do.call("c",do.call("c", output.list))


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
