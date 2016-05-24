##' Get all office types tracked
##' 
##' This function is a wrapper for the Office.getTypes() method of the PVS API Office class which grabs a list of office types Project Vote Smart keeps track of.
##' @usage Office.getTypes()
##' @return A data frame with rows for each office type and columns with the following variables describing the office type:\cr officeTypes.type*.officeTypeId,\cr officeTypes.type*.officeLevelId,\cr officeTypes.type*.officeBranchId,\cr officeTypes.type*.name.
##' @references http://api.votesmart.org/docs/Office.html
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get all office types tracked
##' \dontrun{officetypes <- Office.getTypes()}
##' \dontrun{officetypes}
##' @export


Office.getTypes <-
function () {
  
  
  
  request <-  "Office.getTypes?"
  inputs  <-  ""
  output  <-  pvsRequest4(request,inputs)
  
  output
  
  
}
