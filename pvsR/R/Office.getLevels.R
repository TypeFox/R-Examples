##' Get a list of levels of government and their IDs
##' 
##' This function is a wrapper for the Office.getLevels() method of the PVS API Office class which grabs a list of the levels of government and their IDs.
##' @usage Office.getLevels()
##' @return A data frame with a row for each level and columns with the following variables describing the level:\cr levels.level*.officeLevelId,\cr levels.level*.name.
##' @references http://api.votesmart.org/docs/Office.html
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of government levels and their IDs
##' \dontrun{levels <- Office.getLevels()}
##' \dontrun{levels}

##' @export


Office.getLevels <-
function () {
  
  
  
  request <-  "Office.getLevels?"
  inputs  <-  ""
  output  <-  pvsRequest4(request,inputs)
  
  output
  
  
}
