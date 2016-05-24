#' Delete a campaign from the server
#' @param campaign_urn campaign uuid
#' @param ... other arguments passed to ohmage.
#' @export
oh.campaign.delete <- function(campaign_urn, ...){
	xhr <- oh.call("/campaign/delete",campaign_urn=campaign_urn, ...);	
	message("Campaign deleted!");
}
