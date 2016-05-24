#' Delete users from the server
#' @param user_list comma seperated list of users
#' @param ... other arguments passed on to ohmage
#' @export
oh.user.delete <- function(user_list, ...){
	xhr <- oh.call("/user/delete",user_list=user_list, ...);	
	message("User deleted!");
}
