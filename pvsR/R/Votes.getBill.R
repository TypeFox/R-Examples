##' Get general information on a bill
##' 
##' This function is a wrapper for the Votes.getBill() method of the PVS API Votes class which grabs the general information on a bill. The function sends a request with this method to the PVS API for all bill IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Votes.getBill(billId,...)
##' @param billId a character string or list of character strings with the bill ID(s) (see references for details)
##' @param ... further arguments that are passed on to internal functions. Currently the argument separate can be defined: separate is a vector of character strings defining subnodes that should be returned separately (e.g., "sponsors"). 
##' @return If separate is not specified, a data frame with a row for each bill and columns with variables describing the bill. If separate is specified, a list containing several data frames, one for each subnode mentioned in separate and additionally one data frame containing all remaining nodes not mentioned in separate. The returned data frame contains a row for each bill and columns with the following variables describing the bill:\cr bill.billnumber,\cr bill.parentbill,\cr bill.title,\cr bill.officialtitle,\cr bill.dateintroduced,\cr bill.type,\cr bill.categories.category*.categoryId,\cr bill.categories.category*.name,\cr bill.billtextLink,\cr bill.sponsors.sponsor*.candidateId,\cr bill.sponsors.sponsor*.name,\cr bill.sponsors.sponsor*.type,\cr bill.committeeSponsors.committeeSponsor*.committeeId,\cr bill.committeeSponsors.committeeSponsor*.name,\cr bill.actions.action*.actionId,\cr bill.actions.action*.level,\cr bill.actions.action*.stage,\cr bill.actions.action*.outcome,\cr bill.actions.action*.statusDate,\cr bill.actions.action*.rollNumber,\cr bill.actions.action*.yea,\cr bill.actions.action*.nay,\cr bill.actions.action*.voiceVote,\cr bill.amendments.amendment*.billNumber,\cr bill.amendments.amendment*.actionId,\cr bill.amendments.amendment*.title,\cr bill.amendments.amendment*.statusDate.
##' @references http://api.votesmart.org/docs/Votes.html\cr
##' Use Votes.getByBillNumber(), Votes.getBillsByCategoryYearState(), Votes.getBillsByYearState(), Votes.getBillsByOfficialYearOffice(), Votes.getBillsByOfficialCategoryOffice(), Votes.getByOfficial(), Votes.getBillsBySponsorYear(), Votes.getBillsBySponsorCategory() or Votes.getBillsByStateRecent() to get a list of bill IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get information about certain bills
##' \dontrun{billinfo <- Votes.getBill(list(2819,6427))}
##' \dontrun{billinfo}
##' # let some variables with subnodes be returned separately (here: "sponsors" and "actions")
##' \dontrun{billinfo2 <- Votes.getBill(billId=list(2819,6427,6590),
##' separate=c("sponsors","actions"))}
##' \dontrun{billinfo2}
##' # check the sponsors of the requested bill (argument of separate)... 
##' \dontrun{billinfo2$sponsors}
##' # ... and the usual variables describing the bill (nodes not mentioned in separate)
##' \dontrun{billinfo2$main}
##' @export

Votes.getBill <-
function (billId, ...) {


  
# internal function
Votes.getBill.basic <- function (.billId, ...) {

request <-  "Votes.getBill?"
inputs  <-  paste("&billId=",.billId, sep="")
output  <-  pvsRequest10(request,inputs,...)

# mark output with billId differently depending on if argument "separate" was defined (class of output is different)

if (class(output)=="data.frame") {
  
  output$billId <-.billId
  
} else { # class must be "list", hence separate was defined
  
  for (i in 1:length(output)) output[[i]]["billId"] <- .billId
  
}

output

}


  # Main function
  output.list <- lapply(billId, FUN= function (b) {
    
      Votes.getBill.basic(.billId=b,...)
           
        
    }
  )



# Again, different handling, if seperate was defined. Detecting this is not straightforward, because of possible wrong requests (warning message --> empty data fram returned)
# Hence, check if any of the returned list entries is itself a list


checkclass <- sapply(1:length(output.list), function(j) {
  
  class(output.list[[j]])=="list"
  
})

# if not, proceed as usual:
if (sum(as.numeric(checkclass))==0) {

output.list <- redlist(output.list)
output <- dfList(output.list)

} else { # if there are some lists, proceed differently
  # combine all list entries of the same type in a df and these dfs in a list


nelements <- 1:length(checkclass)
listelements <- nelements[c(checkclass)]
nonlistelements <- nelements[c(!checkclass)]
nodenames <- names(output.list[[listelements[1]]] )

# first process the elements that contain a list
output.list2 <- lapply(nodenames, function(node) {

nodes <- lapply(listelements, function(e) {
  
  output.list[[e]][[node]]
  
})

dfList(nodes)

})

for (i in 1:length(nodenames)) names(output.list2)[i] <- nodenames[i]

# then the remaining, from error messages
output.listerr <- lapply(nonlistelements, function(e){
  
  output.list[[e]]
  
})

output.err <- do.call("rbind",output.listerr)

output.list2$missingData <- output.err

output.list2

} # end else


} # end function

