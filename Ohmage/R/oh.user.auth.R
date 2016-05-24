#' Retrieve hashed passwd for a user
#' @param user username
#' @param password passwd
#' @param ... other argus
#' @return pw hash
#' @export
oh.user.auth <- function(user, password, ...){
	#when using oh.login, argument 'serverurl' is part of the ... ellipse
	xhr <- oh.call("/user/auth", user=user, password=password, ...);
	
	#check response
	if(xhr$result == "success"){
		message("Hashed password retrieved.");
		return(xhr[["hashed_password"]]);
	} 
	return(xhr);
}

