

###########################################################
# extracts used dataset
IRT.data <- function(object, ...) {
    UseMethod("IRT.data")
}
###########################################################
IRT.data.din <- function( object , ... ){
	dat <- object$dat
	attr(dat,"weights") <- object$control$weights
	attr(dat,"group") <- object$control$group
    return(dat)
			}
############################################################			
IRT.data.gdina <- IRT.data.din
IRT.data.gdm <- IRT.data.din
IRT.data.mcdina <- IRT.data.din
IRT.data.slca <- IRT.data.din
#############################################################