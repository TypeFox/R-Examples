## Function Description:
##     This function removes all variables except those which are specified in
##     the given character vector.

rm.all.but <- function(keep=NULL, envir=.GlobalEnv){
	#----[ checking the input ]----#
	## Check the envir attribute
    if(!is.environment(envir)){
    	stop("You should specify an existing environment")
    }
    
    ## ckeck if keep is defined
    if(is.null(keep)){
        stop("The parameter \"keep\" is not defined. It should be a chacter vector with length 1 or more.")
    }
    
	## Check if the keep is a character vector
	if(class(keep)!="character" | typeof(keep)!="character"){
		stop("The value of \"keep\" parameter should be a chacter vector with length 1 or more.")
	}
	
	## check if the length of keep is more than or equal to 1
	if(!length(keep)){
		stop("The value of \"keep\" parameter should be a chacter vector with length 1 or more.")
	}
	
	
	# If the keep has a name that is not in ls(), show error with the location of bad variable name
	if(any(!is.element(keep, ls(envir = as.environment(envir))))){
        bad_var_location <- which(is.element(keep, ls(envir = as.environment(envir)))==FALSE)
        stop(paste("All the items in the keep should be a real existing variable.\nThe item number", bad_var_location, "is not among variables of selected environment!\n", sep=" "))
    }
	

	#----[ processing ]----#
	remove(list = ls(envir)[!(ls(envir) %in% keep)], envir = envir)
}
