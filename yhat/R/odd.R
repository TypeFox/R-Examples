odd <-
function(val) {
########################################################################################
##DESCRIPTION
##Returns true if value is odd; false if true

##REQUIRED ARGUMENTS
##val						Value to check
		
##Returns true if val is odd
	if (((as.integer(val/2))*2)!=val) return(TRUE)	
	return (FALSE)
}

