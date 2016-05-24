##' Get offices tracked according to type
##' 
##' This function is a wrapper for the Office.getOfficesByType() method of the PVS API Office class which grabs a list of offices Project Vote Smart keeps track of according to their type. The function sends a request with this method to the PVS API for all office type IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Office.getOfficesByType(officeTypeId)
##' @param officeTypeId a character string or list of character strings with the office type ID(s) (see references for details)
##' @return A data frame with a row for each office and columns with the following variables describing the office:\cr offices.office*.officeId,\cr offices.office*.officeTypeId,\cr offices.office*.officeLevelId,\cr offices.office*.officeBranchId,\cr offices.office*.name,\cr offices.office*.title,\cr offices.office*.shortTitle.
##' @references http://api.votesmart.org/docs/Office.html\cr
##' See http://api.votesmart.org/docs/semi-static.html or use Office.getTypes() or Office.getOfficesByLevel() to get a list of office types ID(s). 
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get offices tracked for certain office types
##' \dontrun{offices <- Office.getOfficesByType(list("S","K","L"))}
##' \dontrun{head(offices)}
##' @export


Office.getOfficesByType <-
function (officeTypeId) {
  
  
  # internal function
  Office.getOfficesByType.basic <- function (.officeTypeId) {
    
    request <-  "Office.getOfficesByType?"
    inputs  <-  paste("&officeTypeId=",.officeTypeId,sep="")
    output  <-  pvsRequest4(request,inputs)
    
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(officeTypeId, FUN= function (s) {
    Office.getOfficesByType.basic(.officeTypeId=s)
  }
                        )
  
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
