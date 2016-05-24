#' Authenticate with an ohmage server
#' @param user ohmage username
#' @param password ohmage passwd
#' @param serverurl url to the ohmage server
#' @param ... extra parameters for oh.call
#' @importFrom RJSONIO fromJSON
#' @importFrom RCurl postForm getCurlHandle dynCurlReader parseHTTPHeader fileUpload
#' @import methods
#' @export
#' @examples library(Ohmage)
#' \dontrun{
#' #authentication works like a cookie.
#' #oh.login("ohmage.admin", "ohmage.passwd", "https://example.com/app")
#' 
#' #list campaigns you are in
#' #oh.campaign.read()
#' 
#' #read some data
#' #oh.survey_response.read("urn:ohmage:campaign:mycampaign");
#' }
oh.login <- function(user, password, serverurl, ...){
	if(!is.null(getOption("TOKEN"))){
		stop("Already logged in. Please logout first.");
	}
	
	#create a new curl handler for recycling sessions
	options(OhCurlHandle = getCurlHandle());
	options(OhCurlReader = dynCurlReader(getOption("OhCurlHandle"), binary = TRUE));
	
	#try to login
	message("Trying to login to ", serverurl, " with username: ",user);
	mytoken <- oh.auth_token(user, password, serverurl, ...);

	options(SERVERURL = serverurl);
	options(TOKEN = mytoken);
	options(ohmage_username = user);
	
	return(mytoken);
}