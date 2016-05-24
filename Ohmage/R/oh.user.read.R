#' Read info about a user
#' @param user_list username
#' @param ... other args passed to ohmage
#' @export
oh.user.read <- function(user_list, ...){
	xhr <- oh.call("/user/read ", user_list=user_list, ...);		
	return(xhr);	
}
