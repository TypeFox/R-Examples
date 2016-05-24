##' Get a candidate's main biographical information
##' 
##' This function is a wrapper for the CandidateBio.getBio() method of the PVS API CandidateBio class which grabs the main biographical information for each candidate. The function sends a request with this method to the PVS API for all candidate IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage CandidateBio.getBio(candidateId)
##' @param candidateId a character string or list of character strings with the candidate ID(s) (see references for details)
##' @return A data frame with a row for each candidate and columns with the following variables describing the candidate:\cr bio.candidate.crpId (OpenSecrets ID),\cr bio.candidate.firstName,\cr bio.candidate.nickName,\cr bio.candidate.middleName,\cr bio.candidate.lastName,\cr bio.candidate.suffix,\cr bio.candidate.birthDate,\cr bio.candidate.birthPlace,\cr bio.candidate.pronunciation,\cr bio.candidate.gender,\cr bio.candidate.family,\cr bio.candidate.photo,\cr bio.candidate.homeCity,\cr bio.candidate.homeState,\cr bio.candidate.education,\cr bio.candidate.profession,\cr bio.candidate.political,\cr bio.candidate.religion,\cr bio.candidate.congMembership,\cr bio.candidate.orgMembership,\cr bio.candidate.specialMsg,\cr bio.office.parties,\cr bio.office.title,\cr bio.office.shortTitle,\cr bio.office.name,\cr bio.office.type,\cr bio.office.status,\cr bio.office.firstElect,\cr bio.office.lastElect,\cr bio.office.nextElect,\cr bio.office.termStart,\cr bio.office.termEnd,\cr bio.office.district,\cr bio.office.districtId,\cr bio.office.stateId,\cr bio.office.committee*.committeeId,\cr bio.office.committee*.committeeName,\cr bio.election*.office,\cr bio.election*.officeId,\cr bio.election*.officeType,\cr bio.election*.parties,\cr bio.election*.district,\cr bio.election*.districtId,\cr bio.election*.status,\cr bio.election*.ballotName.
##' @references http://api.votesmart.org/docs/CandidateBio.html\cr 
##' Use Candidates.getByOfficeState(), Candidates.getByOfficeTypeState(), Candidates.getByLastname(), Candidates.getByLevenshtein(), Candidates.getByElection(), Candidates.getByDistrict() or Candidates.getByZip() to get a list of candidate IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get main biographical data on Barack Obama and Mitt Romney
##' \dontrun{bio <- CandidateBio.getBio(list(9490,21942))}
##' \dontrun{bio}
##' @export



CandidateBio.getBio <-
function (candidateId) {

  
# internal function (here a special case because several nodes with different number of children are involved, as well as some nodes which only exist for some office types)
CandidateBio.getBio.basic <- function (.candidateId) {
  
request <-  "CandidateBio.getBio?"
inputs  <-  paste("&candidateId=",.candidateId, sep="")
output  <-  pvsRequest3(request,inputs)
output  <-  lapply(output, FUN=function(i){ 
  if (length(grep(pattern="candidateId",x=names(i), value=TRUE))>0) {i}
  else {data.frame(i,"candidateId" =.candidateId)}})
output

}


  # Main function (here also a special case, see above. the output is a list containing dataframes of each main node. each dataframe is processed separately an then merged to one.)
  output.list <- lapply(candidateId, FUN= function (b) {
    
      CandidateBio.getBio.basic(.candidateId=b)
           
        
    }
  )

 nodes <- lapply(output.list,FUN=function(x){names(x)})
 nodes <- unique(unlist(nodes))
 
output.list2 <- lapply(nodes,FUN=function(y) {lapply(output.list, FUN=function(x) {
  df. <-  x[[y]]
  candId <- df.$candidateId
  df.$candidateId <- NULL
  if (!is.null(df.)){
  names(df.) <- paste(y,names(df.),sep=".")
  df. <- data.frame(df.,"candidateId"=candId)
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
  
  candId <- x[,"candidateId"]
  
  x[,"candidateId"] <- NULL
  
  x <- data.frame(cbind(x,matrix(NA, ncol=max.cols-ncol(x)-1, nrow = 1, )),"candidateId"=candId,row.names=NULL)
  
  names(x) <- names(ol[[coln]])
  x
    
} else {x}

}


})

output <- do.call("rbind",ol2)
output

}
                )

max.id <- which.is.max(sapply(output.list3,function(x){dim(na.omit(x["candidateId"]))[1]}));

max.ids <- output.list3[[max.id]]["candidateId"]

output.list4 <- lapply(output.list3, FUN= function(x) {
  x["candidateId"] <- max.ids  
  x
  })

Reduce(function(x,y) {merge(x,y)}, output.list4)

}
