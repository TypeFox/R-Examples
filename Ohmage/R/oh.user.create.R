#' Create a new user	
#' @param username username
#' @param password password
#' @param admin if the user should be an admin
#' @param enabled if new account should be enabled
#' @param new_account user marked as new_account have to change their passwd
#' @param campaign_creation_privilege if user can create campaigns
#' @param ... other arguments passed on to ohmage
#' @export
oh.user.create <- function(username, password, admin = FALSE, enabled = TRUE,  new_account = FALSE, campaign_creation_privilege=FALSE, ...){
	
	admin <- ifelse(admin, "true", "false");
	enabled <- ifelse(enabled, "true", "false");
	new_account <- ifelse(new_account, "true", "false");
	campaign_creation_privilege <- ifelse(campaign_creation_privilege, "true", "false");
	
	xhr <- oh.call("/user/create ", username=username, password=password, admin=admin, enabled=enabled, new_account=new_account, campaign_creation_privilege=campaign_creation_privilege, ...);		
	if(xhr$result == "success") {
		message("User " ,username, " created!");
	} else {
		return(xhr);
	}
}
