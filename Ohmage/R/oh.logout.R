#' Logout from an ohmage session
#' @export
oh.logout <- function(){
	if(is.null(getOption("TOKEN"))){
		stop("Not logged in.");
	}
	options(TOKEN = NULL);
	options(SERVERURL = NULL);
	options(ohmage_username = NULL)
	message("Logged out.\n")
}

