##' Get a list of members of a committee
##' 
##' This function is a wrapper for the Committee.getCommitteeMembers() method of the PVS API Committee class which returns a list of members of a committee. The function sends a request with this method to the PVS API for all committee IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Committee.getCommitteeMembers(committeeId)
##' @param committeeId a character string or list of character strings with the committee ID(s) (see references for details)
##' @return A data frame with a row for each committee member and columns with the following variables describing the committee member:\cr committeeMembers.committee.committeeId,\cr committeeMembers.committee.parentId,\cr committeeMembers.committee.name,\cr committeeMembers.member*.candidateId,\cr committeeMembers.member*.title,\cr committeeMembers.member*.firstName,\cr committeeMembers.member*.middleName,\cr committeeMembers.member*.lastName,\cr committeeMembers.member*.suffix,\cr committeeMembers.member*.party,\cr committeeMembers.member*.position.
##' @references http://api.votesmart.org/docs/Committee.html\cr
##' Use CandidateBio.getBio(), Committee.getCommitteesByTypeState() or Votes.getBill() to get committee ID(s).
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get a list of members of certain committees
##' \dontrun{comember <- Committee.getCommitteeMembers(1)}
##' \dontrun{comember}
##' @export


Committee.getCommitteeMembers <-
function (committeeId) {

  
#internal function:
Committee.getCommitteeMembers.basic <- function (.committeeId) {
  
request <-  "Committee.getCommitteeMembers?"
inputs  <-  paste("&committeeId=",.committeeId, sep="")

pvs.url <- paste("http://api.votesmart.org/",request,"key=",get('pvs.key',envir=.GlobalEnv),inputs,sep="") #generate url for request

doc <- xmlTreeParse(pvs.url)
a <- xmlRoot(doc)

if (length(a)==1 && names(a[1])=="errorMessage") {
  output.items.df <- data.frame("committeeId"=.committeeId)
} else {

items <- getNodeSet(a, path="//member")


output.items <- lapply(items, function(x) data.frame(t(xmlSApply(x, xmlValue))))


output.items.df <- do.call("rbind", output.items)
output.base.df <- data.frame(t(xmlSApply(a[[2]], xmlValue)))
output.df <- merge(output.base.df, output.items.df)


}
output.df

}



  # Main function 
  output.list <- lapply(committeeId, FUN= function (b) {
    
      Committee.getCommitteeMembers.basic(.committeeId=b)
           
        
    }
  )




output <- do.call("rbind", output.list)

output

}
