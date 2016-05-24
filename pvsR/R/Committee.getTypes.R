##' Get the committee types (house, senate, joint, etc.)
##' 
##' This function is a wrapper for the Committee.getTypes() method of the PVS API Committee class which returns the existing committee types.
##' @usage Committee.getTypes()
##' @return A data frame with a row for each committee type and columns with the following variables describing the committee type:\cr committeeTypes.type*.committeeTypeId,\cr committeeTypes.type*.name.
##' @references http://api.votesmart.org/docs/Committee.html
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get the committee types
##' \dontrun{comtypes <- Committee.getTypes()}
##' \dontrun{comtypes}
##' @export



Committee.getTypes <-
function() {
  
  request <-  "Committee.getTypes?"
  inputs  <-  ""
  url.base <- "http://api.votesmart.org/"

  pvs.url <- paste(url.base,request,"key=",get('pvs.key',envir=.GlobalEnv),inputs,sep="") #generate url for request
  output <- t(xmlSApply(removeChildren(xmlRoot(xmlTreeParse(pvs.url,useInternalNodes=TRUE)),kids=1), function(x) xmlSApply(x, xmlValue)))
  data.frame(output, row.names=NULL)
  data.frame(output, row.names=NULL)
  
}
