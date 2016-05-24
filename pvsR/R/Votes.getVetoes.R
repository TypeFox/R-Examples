##' Get a list of vetoes according to candidate
##' 
##' This function is a wrapper for the Votes.getVetoes() method of the PVS API Votes class which returns a list of vetoes according to candidate. The function sends a request with this method to the PVS API for all candidate IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Votes.getVetoes(candidateId)
##' @param candidateId a character string or list of character strings with the candidate ID(s) (see references for details)
##' @return A data frame with a row for each veto and columns with the following variables describing the veto:\cr bills.bill*.vetoId,\cr bills.bill*.statusDate,\cr bills.bill*.billId,\cr bills.bill*.billNumber,\cr bills.bill*.billTitle,\cr bills.bill*.vetoCode,\cr bills.bill*.vetoType,\cr bills.bill*.billSummary,\cr bills.bill*.billLink,\cr bills.bill*.vetoLetterLink.
##' @references http://api.votesmart.org/docs/Votes.html\cr
##' Use Candidates.getByOfficeState(), Candidates.getByOfficeTypeState(), Candidates.getByLastname(), Candidates.getByLevenshtein(), Candidates.getByElection(), Candidates.getByDistrict() or Candidates.getByZip() to get a list of candidate IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get vetoes by Barack Obama
##' \dontrun{vetoes <- Votes.getVetoes(9490)}
##' \dontrun{vetoes}
#' @export


Votes.getVetoes <-
function (candidateId) {
  
  
# internal function

Votes.getVetoes.basic <- function (.candidateId) {
  
request <-  "Votes.getVetoes?"
inputs  <-  paste("&candidateId=",.candidateId, sep="")
output  <-  pvsRequest(request,inputs)
output$candidateId  <- .candidateId  
output

}



# Main function  
  output.list <- lapply(candidateId, FUN= function (b) {
    
      Votes.getVetoes.basic(.candidateId=b)
           
        
    }
  )

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
