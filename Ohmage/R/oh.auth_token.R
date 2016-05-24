oh.auth_token <- function(user, password, ...){
	#when using oh.login, argument 'serverurl' is part of the ... ellipse
	xhr <- oh.call("/user/auth_token", user=user, password=password, ...);
	message("Login successful!");
	return(xhr$token);
}


