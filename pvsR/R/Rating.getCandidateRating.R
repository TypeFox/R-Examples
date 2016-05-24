##' Get a candidate's rating by special interest groups
##' 
##' This function is a wrapper for the Rating.getCandidateRating() method of the PVS API Rating class which grabs a candidate's rating by special interest groups (SIG). The function sends a request with this method to the PVS API for all candidate and SIG IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage Rating.getCandidateRating(candidateId, sigId=NULL)
##' @param candidateId a character string or list of character strings with the candidateId(s) (see references for details)
##' @param sigId (optional) a character string or list of character strings with the special interest group's ID(s) (see references for details)
##' @return A data frame with a row for each rating of a candidate and columns with the following variables describing the candidate:\cr candidateRating.candidate.title,\cr candidateRating.candidate.firstName,\cr candidateRating.candidate.middleName,\cr candidateRating.candidate.lastName,\cr candidateRating.candidate.suffix,\cr candidateRating.candidate.office,\cr candidateRating.rating*.sigId,\cr candidateRating.rating*.ratingId,\cr candidateRating.rating*.categories.category*.categoryId,\cr candidateRating.rating*.categories.category*.name,\cr candidateRating.rating*.timeSpan,\cr candidateRating.rating*.rating,\cr candidateRating.rating*.ratingName,\cr candidateRating.rating*.ratingText.
##' @references http://api.votesmart.org/docs/Rating.html\cr
##' Use Candidates.getByOfficeState(), Candidates.getByOfficeTypeState(), Candidates.getByLastname(), Candidates.getByLevenshtein(), Candidates.getByElection(), Candidates.getByDistrict() or Candidates.getByZip() to get a list of candidate IDs.\cr
##' Use Rating.getSigList() to get a list of special interest group's IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get ratings by candidate and special interest group
##' \dontrun{rating <- Rating.getCandidateRating(candidateId="9490")}
##' \dontrun{head(rating)}
##' @export






Rating.getCandidateRating <-
function (candidateId, sigId=NULL) {
  
  if (length(sigId)==0) {
    
    # internal function
    Rating.getCandidateRating.basic1 <- function (.candidateId) {
      
      request <-  "Rating.getCandidateRating?"
      inputs  <-  paste("&candidateId=",.candidateId,sep="")
      output  <-  pvsRequest6.1b(request,inputs)
      output$candidateId <- .candidateId
      output
      
    }
    
    
    # Main function  
    
    output.list <- lapply(candidateId, FUN= function (y) {
      
        Rating.getCandidateRating.basic1(.candidateId=y)
      
      
    }
    )  
    
    
    
    
    
  } else {
  
  # internal function
  Rating.getCandidateRating.basic2 <- function (.candidateId, .sigId) {
    
    request <-  "Rating.getCandidateRating?"
    inputs  <-  paste("&candidateId=",.candidateId,"&sigId=",.sigId,sep="")
    output  <-  pvsRequest6.1b(request,inputs)
    output$candidateId <- .candidateId
    output
    
  }
  
  
  # Main function  
  
  output.list <- lapply(candidateId, FUN= function (y) {
    lapply(sigId, FUN= function (s) {
      Rating.getCandidateRating.basic2(.candidateId=y, .sigId=s)
    }
           )
  }
                        )
  
  
  
  }
  
  output.list <- redlist(output.list)
  
  
  output <- dfList(output.list)
  
  
  output
  
  
  
}
