## Function Description:
##     This function gets a list of variables as a character vector and save
##     each variable name in a separate file, so that they can be loaded
##     separately.

save.var <- function(varlist=ls(envir = as.environment(.GlobalEnv)), path=getwd(), newdir=TRUE, newdirtag=NULL, envir=.GlobalEnv){
    #----[ checking input ]----#
    ## Check the envir attribute
    if(!is.environment(envir)){
    	stop("You should specify an existing environment")
    }
    
    ## check the varlist attribute
    # Check if the varlist is empty, show error and mention default
    if(!length(varlist)){
       stop("The varlist parameter should contain an existing variable.")
    }
    # If the varlist has a name that is not in ls(), show error with the location of bad variable name
    if(sum(is.element(varlist, ls(envir = as.environment(envir))))!=length(varlist)){
        bad_var_location <- which(is.element(varlist, ls())==FALSE)
        stop(paste("All the items in the varlist should be a real existing variable.\nThe item number", bad_var_location, "is not among variables of selected environment!\n", sep=" "))
    }
    
    ## check the path attribute
    # The path should be a vector of size one and the folder should exist
    if(length(path)!=1){
    	stop("The path attribute should be a vector with length 1.")
    }
    # checking if the path exists
    if (!file.exists(path)){
    	stop("The directory path should fully exist. Please use a correct existing path.")
    }
    
    ## check the newdir attribute
    if(newdir!=TRUE & newdir!=FALSE){
    	stop("The newdir attribute only accepts TRUE or FALSE as input.")
    }
    
    ## check newdirtag attribute
    if(length(newdirtag)){
    	# The newdirtag should be a vector of size one
    	if(length(newdirtag)!=1){
    		stop("The newdirtag attribute should be a vector with length 1.")
    	}
    	# If the newdirtag was empty or contained more than 240 character
        if(nchar(newdirtag)>240){
    	    stop("The acceptable number of characters for newdirtag attribute is 1 up to 240.")
        }
    }
    
    

    #----[ create new directory ]----#
    if(newdir==TRUE){
    	# Getting time
    	time.now <- format(Sys.time(), "%Y%m%d-%H%M%S")
    	if(length(newdirtag)){
    		# Create new folder name with tag
    		final.dir = paste(time.now, "-", newdirtag, sep='')
    	}else{
    		# Create new folder name without tag
    		final.dir = time.now
    	}
    	dir.create(file.path(path, final.dir))
    }

    #----[ hold current folder ]----#
    # save current folder in a variable so that we set it back after the process is finished
    user.wd <- getwd()


    #----[ save variables ]----#
    # Switch to newly created directory
    setwd(file.path(path, final.dir))
    # Write variables one by one in separate files.
    for(i in varlist){
        save(file=paste(i, ".RData"), list=i, envir = envir)
    }

    #----[ return to original folder ]----#
    # Switch back to user's working directory
    setwd(user.wd)
}

