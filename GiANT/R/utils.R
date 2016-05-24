##########################################################################
#getNeededParameters
#
# f_func:			A list of all parameters needed for a function as
#					given by 'formals'.
# given_parameters:	A list of parameters given to the main function.
#
# The function checks if there are more parameters needed for a function
# than given. If these parameters are predefined the function adds the
# missing parameters to the list of paramters which is returned.
##########################################################################
getNeededParameters <- function(f_func, given_parameters){
	f_nms <- names(f_func)
	given_parameters <- given_parameters[!is.na(names(given_parameters))]

	ind <- sapply(f_nms, function(x){
			if(length(names(given_parameters)) == 0){
				if(is.null(f_func[[x]])){
					return(TRUE)
				}else if(f_func[[x]] == ""){
					return(FALSE)
				}else{
					return(TRUE)
				}
			}else{
				if(x %in% names(given_parameters)){
					return(FALSE)
				}else if(is.null(f_func[[x]])){
					return(TRUE)
				}else if(f_func[[x]] == ""){
					return(FALSE)
				}else{
					return(TRUE)
				}
			}
		})

	if(length(f_func[ind]) == 0){
		return(given_parameters)
	}
	if(length(names(given_parameters)) == 0){
		return(f_func[ind])
	}else{
		if(is.na(names(given_parameters)[1])){
			return(f_func[ind])
		}else{
			return(c(f_func[ind], given_parameters))
		}
	}
}