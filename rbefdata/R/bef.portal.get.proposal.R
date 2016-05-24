#' Fetch primary data in CSV format from a BEFdata portal paperproposal
#'
#' This function fetches data associated to a paperproposal on a BEFdata
#' portal.  By default it will fetch all associated files in CSV format next to
#' the informations available like the title, the download URL and the id of
#' the datasets.
#'
#' You need to provide the function with a proposal id which you can find in
#' the URL of the proposal and your user credentials. You can find the credentials
#' inside of your profile page on the BEFdata portal. The credentials ensure you
#' have the rights to download the data.
#'
#' The function returns a list object which you can store to a variable as
#' shown in the examples below.
#'
#' @param id This is the ID of a paper proposal. You can download all datasets in
#'        one turn given the ID.
#' @param curl If using in a loop, call getCurlHandle() first and pass
#'        the returned value in here (avoids unnecessary footprint)
#' @param \dots This are other arguments passed to \code{\link[RCurl]{getURLContent}}
#'
#' @return The function returns a list of raw data attached to a proposal.
#' 	   An error is thrown when the proposal is not found or you don't have
#'	   the access rights for it.
#'
#' @examples \dontrun{
#'	  prop1 = bef.portal.get.datasets.for_proposal(proposal = 8)
#'  	}
#' @import RCurl
#' @export bef.portal.get.datasets.for_proposal bef.get.datasets.for_proposal
#' @aliases bef.get.datasets.for_proposal

bef.portal.get.datasets.for_proposal <- bef.get.datasets.for_proposal <- function(id, curl = getCurlHandle(), ...) {
  paperproposal_url = paperproposal_url(proposal_id = id)
  proposal_raw_csv = getURLContent(paperproposal_url, curl = curl, ...)
  if (getCurlInfo(curl)$response.code != 200) {
    stop("Proposal not found or not accessible. Please check your credentials and make sure you have access right for it.")
  }
  proposal_data = read.csv(text = proposal_raw_csv)
  datasets = lapply(proposal_data$ID, function(x) bef.portal.get.dataset(id = x))
  return(datasets)
}
