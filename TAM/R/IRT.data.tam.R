
##############################################
# IRT.data functions
IRT.data.tam.mml <- function( object , ... ){
	dat <- object$resp
	attr(dat,"weights") <- object$pweights
	attr(dat,"group") <- object$group
    return(dat)
			}
IRT.data.tam.mml.3pl <- IRT.data.tam.mml			
IRT.data.tamaan <- IRT.data.tam.mml
#################################################