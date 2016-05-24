##' Get a list of branches of government and their IDs
##' 
##' This function is a wrapper for the Office.getBranches() method of the PVS API Office class which grabs a list of branches of government and their IDs.
##' @usage Office.getBranches()
##' @return A data frame with a row for each branch and columns with the following variables describing the branch:\cr branches.branch*.officeBranchId,\cr branches.branch*.name.
##' @references http://api.votesmart.org/docs/Office.html
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of government branches and their IDs
##' \dontrun{branches <- Office.getBranches()}
##' \dontrun{branches}

##' @export


Office.getBranches <-
function () {
  
  
  
  request <-  "Office.getBranches?"
  inputs  <-  ""
  output  <-  pvsRequest4(request,inputs)
  
  output
  
  
}
