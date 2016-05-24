##' Get a list of bills according to office (optional), candidate and year
##' 
##' This function is a wrapper for the Votes.getBillsByOfficialYearOffice() method of the PVS API Votes class which grabs a list of bills that fit the candidate and year. The function sends a request with this method to the PVS API for all years, candidate and office IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Votes.getBillsByOfficialYearOffice(year, candidateId, officeId=NULL)
##' @param year a character string or list of character strings with the year ID(s)
##' @param candidateId a character string or list of character strings with the candidate ID(s) (see references for details)
##' @param officeId (optional) a character string or list of character strings with the office ID(s) (default: all) (see references for details)
##' @return A data frame with a row for each bill and columns with the following variables describing the bill:\cr bills.bill*.billId,\cr bills.bill*.billNumber,\cr bills.bill*.title,\cr bills.bill*.type.
##' @references http://api.votesmart.org/docs/Votes.html\cr
##' Use Candidates.getByOfficeState(), Candidates.getByOfficeTypeState(), Candidates.getByLastname(), Candidates.getByLevenshtein(), Candidates.getByElection(), Candidates.getByDistrict() or Candidates.getByZip() to get a list of candidate IDs.\cr
##' See http://api.votesmart.org/docs/semi-static.html for a list of office IDs or use Office.getOfficesByType(), Office.getOfficesByLevel(), Office.getOfficesByTypeLevel() or Office.getOfficesByBranchLevel
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of bills for a certain office, candidate and year 
##' \dontrun{bills <- Votes.getBillsByOfficialYearOffice(list(2010,2011,2012),107800,)}
##' \dontrun{bills}
##' @export




Votes.getBillsByOfficialYearOffice <-
function (year, candidateId, officeId=NULL) {
  
    if (length(officeId)==0) {
      # internal function
Votes.getBillsByOfficialYearOffice.basic1 <- function (.year, .candidateId) {
  
request <-  "Votes.getBillsByOfficialYearOffice?"
inputs  <-  paste("&year=",.year,"&candidateId=",.candidateId,sep="")
output  <-  pvsRequest(request,inputs)
output$year <- .year
output$candidateId <- .candidateId
output

}
  

# Main function  

  output.list <- lapply(year, FUN= function (y) {
    lapply(candidateId, FUN= function (s) {
      Votes.getBillsByOfficialYearOffice.basic1(.year=y, .candidateId=s)
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
      
      
    } else {
      
# internal function
Votes.getBillsByOfficialYearOffice.basic2 <- function (.year, .candidateId, .officeId) {
  
request <-  "Votes.getBillsByOfficialYearOffice?"
inputs  <-  paste("&year=",.year,"&candidateId=",.candidateId, "&officeId=", .officeId, sep="")
output  <-  pvsRequest(request,inputs)
output$year <- .year
output$candidateId <- .candidateId
output$officeId.input <- .officeId
output

}
  

# Main function  

  output.list <- lapply(year, FUN= function (y) {
    lapply(candidateId, FUN= function (s) {
      lapply(officeId, FUN= function (c) {
       Votes.getBillsByOfficialYearOffice.basic2(.year=y, .candidateId=s, .officeId=c)
              }
             )
           }
        )
    }
  )
  
  
# reduce lists in list (3fold) to one list
output.list <- do.call("c",do.call("c", output.list))


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

# Avoids, that output is missleading, because officeId is already given in request-output, but also a
# additionally generated (as officeId.input). Problem exists because some request-outputs might be empty
# and therefore only contain one "officeId" whereas the non-empty ones contain two. (see basic function)
output$officeId[c(as.vector(is.na(output$officeId)))] <- output$officeId.input[as.vector(is.na(output$officeId))]
output$officeId.input <- NULL

}

output

}
