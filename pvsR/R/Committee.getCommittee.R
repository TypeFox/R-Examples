##' Get detailed committee (contact) information
##' 
##' This function is a wrapper for the Committee.getCommittee() method of the PVS API Committee class which returns detailed committee data. The function sends a request with this method to the PVS API for all committee IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Committee.getCommittee(committeeId)
##' @param committeeId a character string or list of character strings with the committee ID(s) (see references for details)
##' @return A data frame with a row for each committee and columns with the following variables describing the committee:\cr committee.committeeId,\cr committee.parentId,\cr committee.stateId,\cr committee.committeeTypeId,\cr committee.name,\cr committee.jurisdiction,\cr committee.contact.address1,\cr committee.contact.address2,\cr committee.contact.city,\cr committee.contact.state,\cr committee.contact.zip,\cr committee.contact.phone,\cr committee.contact.fax,\cr committee.contact.email,\cr committee.contact.url,\cr committee.contact.staffContact.
##' @references http://api.votesmart.org/docs/Committee.html\cr 
##' Use CandidateBio.getBio(), Committee.getCommitteesByTypeState() or Votes.getBill() to get committee ID(s).
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get information about a certain committee
##' \dontrun{committee <- Committee.getCommittee(list(1,2))}
##' \dontrun{committee}
##' @export


Committee.getCommittee <-
function (committeeId) {

  
# internal function (here a special case because several nodes with different number of children are involved, as well as some nodes which only exist for some office types)
Committee.getCommittee.basic <- function (.committeeId) {
  
request <-  "Committee.getCommittee?"
inputs  <-  paste("&committeeId=",.committeeId, sep="")
output  <-  pvsRequest3(request,inputs)
output  <-  lapply(output, FUN=function(i){ 
  if (length(grep(pattern="committeeId",x=names(i), value=TRUE))>0) {i}
  else {data.frame(i,"committeeId" =.committeeId)}})
output

}


  # Main function (here also a special case, see above. the output is a list containing dataframes of each main node. each dataframe is processed separately an then merged to one.)
  output.list <- lapply(committeeId, FUN= function (b) {
    
      Committee.getCommittee.basic(.committeeId=b)
           
        
    }
  )

 nodes <- lapply(output.list,FUN=function(x){names(x)})
 nodes <- unique(unlist(nodes))
 
output.list2 <- lapply(nodes,FUN=function(y) {lapply(output.list, FUN=function(x) {
  df. <-  x[[y]]
  candId <- df.$committeeId
  df.$committeeId <- NULL
  if (!is.null(df.)){
  names(df.) <- paste(y,names(df.),sep=".")
  df. <- data.frame(df.,"committeeId"=candId)
  df.
  }
  })})

output.list3 <- lapply(output.list2, FUN=function(ol){

# which list entry has the most columns, how many are these?
coln <- which.is.max(sapply(ol,function(x){if (!is.null(x)){ncol(x)} else {0}}));
max.cols <- max(sapply(ol,function(x){if (!is.null(x)){ncol(x)} else {0}}));

# give all list entries (dfs in list) the same number of columns and the same names
ol2 <- lapply(ol, function(x){

if (is.null(x)) {NA}

else {
  
if (ncol(x) < max.cols) {
  
  candId <- x[,"committeeId"]
  
  x[,"committeeId"] <- NULL
  
  x <- data.frame(cbind(x,matrix(NA, ncol=max.cols-ncol(x)-1, nrow = 1, )),"committeeId"=candId,row.names=NULL)
  
  names(x) <- names(ol[[coln]])
  x
    
} else {x}

}


})

output <- do.call("rbind",ol2)
output

}
                )

max.id <- which.is.max(sapply(output.list3,function(x){dim(na.omit(x["committeeId"]))[1]}));

max.ids <- output.list3[[max.id]]["committeeId"]

output.list4 <- lapply(output.list3, FUN= function(x) {
  x["committeeId"] <- max.ids  
  x
  })

output <- Reduce(function(x,y) {merge(x,y)}, output.list4)

names(output) <- sub(pattern=".text", replacement="", x=names(output), fixed=TRUE)

output

}
