ENmiscEnv <- function() {
#	pos <-  match("ENmiscEnv", search()) 
#	if (is.na(pos)) {
#		ENmiscEnv <- list()
#		attach(ENmiscEnv, pos = length(search()) - 1) 
#		rm (ENmiscEnv) 
#		pos <- match("ENmiscEnv", search()) 
#	} 
#	return(pos.to.env(pos)) 
	return(ENmiscEnvironment)
}


putENmisc <- function(x, value) assign(x, value, envir = ENmiscEnv())

getENmisc <- function(x, mode="any") get(x, envir = ENmiscEnv(), mode = mode, inherits = FALSE)

ENmiscEnvironment <- new.env()

#####################################

# putENmisc("putRENmisc",putENmisc)

# putENmisc("getENmisc",getENmisc)

# putENmisc("ENmiscEnv",ENmiscEnv)

# rm(getENmisc,putENmisc,ENmiscEnv)

