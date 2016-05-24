##' Get district IDs according to the office and state
##' 
##' This function is a wrapper for the District.getByOfficeState() method of the PVS API District class which grabs district IDs according to the office and state. The function sends a request with this method to the PVS API for all district names and office IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage District.getByOfficeState(officeId, stateId, districtName=NULL, all=FALSE)
##' @param officeId a character string or list of character strings with the office ID(s) (see references for details)
##' @param stateId a character string or list of character strings with the state ID(s) (see references for details)
##' @param districtName (optional) a character string or list of character strings with the district name (default: All)
##' @param all a logical indicator; if TRUE data on all possible combinations of the input variables are returned, if FALSE (default) only the exact combinations of them (see example)
##' @return A data frame with a row for each district and columns with the following variables describing the district:\cr districtList.district*.districtId,\cr districtList.district*.name,\cr districtList.district*.officeId,\cr districtList.district*.stateId.
##' @references http://api.votesmart.org/docs/District.html\cr
##' See http://api.votesmart.org/docs/semi-static.html for a list of office IDs or use Office.getOfficesByType(), Office.getOfficesByLevel(), Office.getOfficesByTypeLevel() or Office.getOfficesByBranchLevel().\cr
##' Use State.getStateIDs() to get a list of state IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get districts for certain office and state IDs 
##' \dontrun{districts <- District.getByOfficeState(officeId=list(8,9),stateId=list("NY","NJ"))}
##' \dontrun{district}
##' # get a data frame of districts according to all state/office/districtName combinations
##' \dontrun{districts <- District.getByOfficeState(officeId=list(8,9),stateId=list("NY","NJ"),
##'  districtName=list(1,2), all=TRUE)}
##' \dontrun{districts}
##' # get a data frame of districts according to the exact state/office/districtName combinations 
##' # (i.e., 8/"NY"/1, 9/"NJ"/2)
##' \dontrun{districts <- District.getByOfficeState(officeId=list(8,9),stateId=list("NY","NJ"),
##'  districtName=list(1,2), all=FALSE)}
##' \dontrun{districts}
##' @export


District.getByOfficeState <-
function (officeId, stateId, districtName=NULL, all=FALSE) {
  
    if (length(districtName)==0) {
      # internal function
District.getByOfficeState.basic1 <- function (.officeId, .stateId) {
  
request <-  "District.getByOfficeState?"
inputs  <-  paste("&officeId=",.officeId,"&stateId=",.stateId,sep="")
output  <-  pvsRequest(request,inputs)
output$officeId <- .officeId
output$stateId <- .stateId
output

}
  

if (all==TRUE) {

# Main function  

  output.list <- lapply(officeId, FUN= function (y) {
    lapply(stateId, FUN= function (s) {
      District.getByOfficeState.basic1(.officeId=y, .stateId=s)
           }
        )
    }
  )
  
  
} else {
  
  # Main function  
  
  reqdf <- data.frame(o=unlist(officeId), s=unlist(stateId))
  
  output.list <- lapply(1:dim(reqdf)[1], FUN= function (l) {
    
    District.getByOfficeState.basic1(.stateId=reqdf[l,"s"], .officeId=reqdf[l,"o"])
    
  })
  
}

output.list <- redlist(output.list)


output <- dfList(output.list)

      
    } else {
      
# internal function
District.getByOfficeState.basic2 <- function (.officeId, .stateId, .districtName) {
  
request <-  "District.getByOfficeState?"
inputs  <-  paste("&officeId=",.officeId,"&stateId=",.stateId, "&districtName=", .districtName, sep="")
output  <-  pvsRequest(request,inputs)
output$officeId <- .officeId
output$stateId <- .stateId
output$districtName.input <- .districtName
output

}
  

if (all==TRUE) {
  

# Main function  

  output.list <- lapply(officeId, FUN= function (y) {
    lapply(stateId, FUN= function (s) {
      lapply(districtName, FUN= function (c) {
       District.getByOfficeState.basic2(.officeId=y, .stateId=s, .districtName=c)
              }
             )
           }
        )
    }
  )
  
  
} else {
  
  # Main function  
  
  reqdf <- data.frame(o=unlist(officeId), s=unlist(stateId), n=unlist(districtName))
  
  output.list <- lapply(1:dim(reqdf)[1], FUN= function (l) {
    
    District.getByOfficeState.basic2(.stateId=reqdf[l,"s"], .officeId=reqdf[l,"o"], .districtName=reqdf[l,"n"])
    
  })
  
}

output.list <- redlist(output.list)

output <- dfList(output.list)

# Avoids, that output is missleading, because districtName is already given in request-output, but also a
# additionally generated (as districtName.input). Problem exists because some request-outputs might be empty
# and therefore only contain one "districtName" whereas the non-empty ones contain two. (see basic function)
output$districtName[c(as.vector(is.na(output$districtName)))] <- output$districtName.input[as.vector(is.na(output$districtName))]
output$districtName.input <- NULL

}
   

output

}
