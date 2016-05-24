##' Get basic data on all offices
##' 
##' This function is essentially a wrapper around Office.getOfficesByLevel().
##' @usage getOffices()
##' @return A data frame with a row for each office and columns with the following variables describing the office:\cr offices.office*.officeId,\cr offices.office*.officeTypeId,\cr offices.office*.officeLevelId,\cr offices.office*.officeBranchId,\cr offices.office*.name,\cr offices.office*.title,\cr offices.office*.shortTitle.
##' @references http://api.votesmart.org/docs/Office.html
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of all offices on all levels
##' \dontrun{offices <- getOffices()}
##' \dontrun{head(offices)}
##' @export

getOffices <- 
  function(){
    
    o.list <- lapply(list("F","S", "L"), FUN=function(i){
      
      Office.getOfficesByLevel(i)
      
    })
    
    dfList(o.list)
    
  }

