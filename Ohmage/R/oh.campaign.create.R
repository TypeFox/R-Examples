#' Create a new campaign
#' @param xml campaign.xml
#' @param running_state if campaign is running or not
#' @param privacy_state private or shared
#' @param class_urn_list classes to add
#' @param description a description
#' @param ... other stuff passed to the server
#' @export
oh.campaign.create <- function(xml, running_state = 'running', privacy_state='shared', class_urn_list='', description="My campaign.", ...){
	xhr <- oh.call("/campaign/create", style="httppost", xml=xml, running_state=running_state, privacy_state=privacy_state, class_urn_list=class_urn_list, description=description, ...);
	message("Campaign created!")
}
