##' Get a list of recent bills according to the state.
##' 
##' This function is a wrapper for the Votes.getBillsByStateRecent() method of the PVS API Votes class which returns a list of recent bills according to the state. The maximum number of bills returned is 100. The function sends a request with this method to the PVS API for all states given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Votes.getBillsByStateRecent(amount, state="NA")
##' @param amount (optional) the amount of bills returned (default: 100, max: 100)
##' @param state (optional) a character string or list of character strings with the state(s) (default: "NA", for national) (see references for details)
##' @return A data frame with a row for each bill and columns with the following variables describing the bill:\cr bills.bill*.billId,\cr bills.bill*.billNumber,\cr bills.bill*.title,\cr bills.bill*.type.
##' @references http://api.votesmart.org/docs/Votes.html\cr
##' Use State.getStateIDs() to get a list of state IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of recent bills according to the state
##' \dontrun{recentbills <- Votes.getBillsByStateRecent(40,list("FL","NY"))}
##' \dontrun{recentbills}

##' @export


Votes.getBillsByStateRecent <-
function (amount=100, state="NA") {
  
    
# internal function
Votes.getBillsByStateRecent.basic <- function (.amount, .state) {
  
request <-  "Votes.getBillsByStateRecent?"
inputs  <-  paste("&amount=",.amount,"&state=",.state,sep="")
output  <-  pvsRequest(request,inputs)
output$state <- .state
output

}
  

# Main function  

  output.list <- lapply(amount, FUN= function (y) {
    lapply(state, FUN= function (s) {
      Votes.getBillsByStateRecent.basic(.amount=y, .state=s)
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
