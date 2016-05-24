ruuid <- function(){
  hex <- c(0:9, letters[1:6]);
  index <- ceiling(runif(32,0,16));
  sym <- hex[index];
  myuuid <- paste(sapply(list(sym[1:8], sym[9:12], sym[13:16], sym[17:20],sym[21:31]), paste, collapse=""), collapse="-");
  return(myuuid);
}


#' Upload a survey response
#' @param campaign_urn campaign id
#' @param user username
#' @param password password
#' @param surveys surveys
#' @param campaign_creation_timestamp ISO 8601 
#' @param ... other arguments passed to ohmage
#' @return data frame with responses
#' @export
oh.survey.upload <- function(campaign_urn, user, password, surveys, campaign_creation_timestamp, ...){
  myuuid <- ruuid();
	xhr <- oh.call("/survey/upload", token=NULL, style="httppost", campaign_urn=campaign_urn, user=user, 
			password=password, surveys=surveys, campaign_creation_timestamp=campaign_creation_timestamp, ...);	
	message("Survey data uploaded for user ", user)	
}
