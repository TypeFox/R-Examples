##' Get a candidate's additional biographical information
##' 
##' This function is a wrapper for the CandidateBio.getAddlBio() method of the PVS API CandidateBio class which grabs the extended biographical information for each candidate that has one. The function sends a request with this method to the PVS API for all candidate-IDs given as a function input, extracts the XML values from the returned XML file(s) and returns them arranged in one data frame.
##' @usage CandidateBio.getAddlBio(candidateId)
##' @param candidateId a character string or list of character strings with the candidate ID(s) (see references for details)
##' @return A data frame with a row for each candidate and columns with the following variables describing the candidate:\cr addlBio.candidate.shortTitle,\cr addlBio.candidate.firstName,\cr addlBio.candidate.nickName,\cr addlBio.candidate.middleName,\cr addlBio.candidate.lastName,\cr addlBio.candidate.suffix,\cr addlBio.additional.item*.name,\cr addlBio.additional.item*.data
##' @references http://api.votesmart.org/docs/CandidateBio.html\cr 
##' Use Candidates.getByOfficeState(), Candidates.getByOfficeTypeState(), Candidates.getByLastname(), Candidates.getByLevenshtein(), Candidates.getByElection(), Candidates.getByDistrict() or Candidates.getByZip() to get a list of candidate IDs.
##' @author Ulrich Matter <ulrich.matter-at-unibas.ch>
##' @examples
##' # First, make sure your personal PVS API key is saved as character string in the pvs.key variable:
##' \dontrun{pvs.key <- "yourkey"}
##' # get additional biographical data on Barack Obama 
##' \dontrun{obama <- CandidateBio.getAddlBio(9490)}
##' \dontrun{obama}
##' # get additional biographical data on Barack Obama and Mitt Romney
##' \dontrun{onr <- CandidateBio.getAddlBio(list(9490,21942))}
##' \dontrun{onr}
##' @export



CandidateBio.getAddlBio <-
function (candidateId) {

  
#internal function:
  CandidateBio.getAddlBio.basic <- function (.candidateId) {
  
request <-  "CandidateBio.getAddlBio?"
inputs  <-  paste("&candidateId=",.candidateId, sep="")
url.base <- "http://api.votesmart.org/"

pvs.url <- paste(url.base,request,"key=",get('pvs.key',envir=.GlobalEnv),inputs,sep="") #generate url for request

doc <- xmlTreeParse(pvs.url)
a <- xmlRoot(doc)

if (length(a)==1 && names(a[1])=="errorMessage") {
  output.items.df <- data.frame("candidateId"=.candidateId)
} else {

items <- getNodeSet(a, path="//item")


output.items <- lapply(items, function(x) data.frame(t(xmlSApply(x, xmlValue))))


output.items <- lapply(output.items, function(x) {
  x. <- data.frame(x["data"])
                  names(x.) <- as.character(x[1,"name"])
                  x.
                  }
       )


output.items.df <- do.call("cbind", output.items)
output.items.df <- data.frame("candidateId"=.candidateId, output.items.df)
}
output.items.df

}



  # Main function (here also a special case, see above. the output is a list containing dataframes of each main node. each dataframe is processed separately an then merged to one.)
  output.list <- lapply(candidateId, FUN= function (b) {
    
      CandidateBio.getAddlBio.basic(.candidateId=b)
           
        
    }
  )



# extract all possible bio vars 

all.bio.vars <- names(do.call("cbind", output.list))

all.bio.vars <- unique(all.bio.vars)


# complete and clean each pm's df with NA's if a var is missing
output.list2 <- lapply(output.list, FUN=function(x) {
  for (i in all.bio.vars ) {
		if ( length(grep(i,names(x), fixed=TRUE))==0 ) {
			x.names.old <- names(x)
			x <- data.frame(x,i=NA)
			names(x) <- c(x.names.old,i)
			
		} else {
			x <- x
		}
	}
	x <- x[,c(all.bio.vars)]
}

)

output <- do.call("rbind", output.list2)

output

}
