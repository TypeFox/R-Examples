##' Get a candidate's detailed biographical information
##' 
##' This function is a wrapper for the CandidateBio.getDetailedBio() method of the PVS API CandidateBio class which grabs the detailed biographical information for each candidate. The function sends a request with this method to the PVS API for all candidate IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage CandidateBio.getDetailedBio(candidateId)
##' @param candidateId a character string or list of character strings with the candidate ID(s) (see references for details)
##' @return A list with several data frames containing the elements of CandidateBio.getBio(), and expands upon:\cr bio.candidate.education.degree,\cr bio.candidate.education.field,\cr bio.candidate.education.school,\cr bio.candidate.education.span,\cr bio.candidate.education.gpa,\cr bio.candidate.education.fullText,\cr bio.candidate.profession.title,\cr bio.candidate.profession.organization,\cr bio.candidate.profession.span,\cr bio.candidate.profession.special,\cr bio.candidate.profession.district,\cr bio.candidate.profession.fullText,\cr bio.candidate.political.title,\cr bio.candidate.political.organization,\cr bio.candidate.political.span,\cr bio.candidate.political.special,\cr bio.candidate.political.district,\cr bio.candidate.political.fullText,\cr bio.candidate.congMembership.title,\cr bio.candidate.congMembership.organization,\cr bio.candidate.congMembership.span,\cr bio.candidate.congMembership.special,\cr bio.candidate.congMembership.district,\cr bio.candidate.congMembership.fullText,\cr bio.candidate.orgMembership.title,\cr bio.candidate.orgMembership.organization,\cr bio.candidate.orgMembership.span,\cr bio.candidate.orgMembership.special,\cr bio.candidate.orgMembership.district,\cr bio.candidate.orgMembership.fullText.
##' @references http://api.votesmart.org/docs/CandidateBio.html\cr 
##' Use Candidates.getByOfficeState(), Candidates.getByOfficeTypeState(), Candidates.getByLastname(), Candidates.getByLevenshtein(), Candidates.getByElection(), Candidates.getByDistrict() or Candidates.getByZip() to get a list of candidate IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get main biographical data on Barack Obama and Mitt Romney
##' \dontrun{bio <- CandidateBio.getDetailedBio(list(9490,21942))}
##' \dontrun{head(bio$profession)}
##' \dontrun{head(bio$orgMembership)}

##' @export



CandidateBio.getDetailedBio <-
function (candidateId) {

  
# internal function (here a special case because several nodes with different number of children are involved, as well as some nodes which only exist for some office types)
CandidateBio.getDetailedBio.basic <- function (.candidateId) {
  
request <-  "CandidateBio.getDetailedBio?"
inputs  <-  paste("&candidateId=",.candidateId, sep="")
output  <-  pvsRequestDetailedBio(request,inputs)
 if (is.data.frame(output)) { 
output <- list(
              candidate=data.frame("candidateId" =.candidateId),
              office=data.frame(name=NA),
              education=data.frame(degree=NA),
              profession=data.frame(title=NA),
              political=data.frame(title=NA),
              congMembership=data.frame(title=NA),
              orgMembership=data.frame(title=NA)
              )
            }
output

}


  # get all output lists of all candidate ids 
  output.list <- lapply(candidateId, FUN= function (b) {
    
      cand.b <- CandidateBio.getDetailedBio.basic(.candidateId=b)
      suppressWarnings(for (j in 1:length(cand.b)) cand.b[[j]]["candidateId"] <- b) # add id to each data frame in list
      cand.b
           
    }
  )

 # combine all data frames of the same category 
 
 output.list2 <- unlist(output.list, recursive=FALSE)
 itemnames <- names(output.list2)
 catnames <- unique(itemnames)

output.list3 <- lapply(catnames, function(i){
  
 cat.i <- data.frame(dfList(output.list2[itemnames %in% i]), row.names=NULL, stringsAsFactors=FALSE)
 for (j in names(cat.i)) cat.i[,j] <- as.character(cat.i[,j])
 return(cat.i)
  
})

names(output.list3) <- catnames


return(output.list3)


}

