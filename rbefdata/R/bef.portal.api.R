#' Basic access to the BEFdata portal XML API.
#'
#' This function gives access to the meta data available for a dataset
#'
#' This function takes a dataset id and returns the metadata of the dataset as list.
#' This is not the eml metadata. It is the metadata read from the XML API of the
#' data portal.
#'
#' @param id Is the ID of the dataset you like to get the information for. You find
#'  	  the ID in the URL of the dataset on the BEFdata portal.
#' @return Returns a list of metadata information for a dataset.
#' @import RCurl
#' @import XML
#' @export bef.portal.api.dataset_info

bef.portal.api.dataset_info <- function(id) {
  dataset_xml_to_list = xmlToList(dataset_url(dataset_id = id, type = "xml"))
  return(dataset_xml_to_list)
}

#' This function returns the information available for a paper proposal.
#'
#' @param id Is the ID of the paperproposal you like to get information for. You can
#'        find the ID in the URL on the paper proposal.
#' @return Returns a list of metadata information for the proposal.
#' @export bef.portal.api.proposal_info

bef.portal.api.proposal_info <- function(id) {
  proposal_xml_to_list = xmlToList(paperproposal_url(proposal_id = id, type = "xml"))
  return(proposal_xml_to_list)
}
