##' Get all the bills an official has voted on by year
##' 
##' This function is a wrapper for the Votes.getByOfficial() method of the PVS API Votes class which dumps all the bills an official has voted on based on the candidateId and year. The function sends a request with this method to the PVS API for all candidate IDs and years given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Votes.getByOfficial(year, candidateId)
##' @param year a character string or list of character strings with the year(s)
##' @param candidateId a character string or list of character strings with the candidate ID(s) (see references for details)
##' @return A data frame with a row for each bill and columns with the following variables describing the bill:\cr bills.bill*.billId,\cr bills.bill*.billNumber,\cr bills.bill*.title,\cr bills.bill.categories.category*.categoryId,\cr bills.bill.categories.category*.name,\cr bills.bill*.officeId,\cr bills.bill*.office,\cr bills.bill*.vote,\cr bills.bill*.actionId,\cr bills.bill*.stage.
##' @references http://api.votesmart.org/docs/Votes.html\cr
##' Use Candidates.getByOfficeState(), Candidates.getByOfficeTypeState(), Candidates.getByLastname(), Candidates.getByLevenshtein(), Candidates.getByElection(), Candidates.getByDistrict() or Candidates.getByZip() to get a list of candidate IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get additional biographical data on Barack Obama 
##' \dontrun{votes <- Votes.getByOfficial(as.list(2010:2012),9490)}
##' \dontrun{votes}

##' @export



Votes.getByOfficial <-
function (year, candidateId) {
  

# internal function
Votes.getByOfficial.basic <- function (.year, .candidateId) {
  
request <-  "Votes.getByOfficial?"
inputs  <-  paste("&year=",.year,"&candidateId=",.candidateId,sep="")
output  <-  pvsRequest(request,inputs)
output$year <-.year
output$candidateId <- .candidateId
output

}


  #Main function
  output.list <- lapply(year, FUN= function (y) {
    lapply(candidateId, FUN= function (s) {
      Votes.getByOfficial.basic(.year=y, .candidateId=s)
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
