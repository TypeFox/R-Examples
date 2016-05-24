#TODO change names, too short #TODO replaced by etaval?
###################################################################################
# Computes the eta value of a signal
#
# @param signal 	The columns of signal correspond to different input components. Must be normalized (zero mean, unit variance)
# @param T			Time interval
#
# @return returns the eta value of the signal in a time interval T time units long.
# @references  \code{\link{leta}}
# @keywords internal
###################################################################################
eta <- function(signal, T){
	res=var(sfaTimediff(signal))
	res=sqrt(res)%*%T/(2*pi)
}


###################################################################################
#' Computes the eta value of a signal (slowness)
#'
#' @param x 	The columns of signal correspond to different input components. Must be normalized (zero mean, unit variance)
#' @param T		Time interval
#'
#' @return returns the eta value of the signal in a time interval T time units long.
# @references  \code{\link{leta}}
#' @export
###################################################################################
etaval <- function(x,T=length(x)){
	return(as.numeric((T/(2*pi))%*%sqrt(mean(diff(x)^2))/(sqrt(mean((x-mean(x))^2)))))
}

###################################################################################
# Compute the eta values of long data signals.
#
# Updates the internal structures and returns
# the eta values of the data signals seen so far.
# The input signal does not need to be normalized.
#
# @param DATA 		The columns of DATA correspond to different input components. Must be normalized (zero mean, unit variance)
# @param t			Time interval
#
# @return returns the eta value of the signal DATA in a time interval T time units long.
#
# @references  \code{\link{eta}}
# @keywords internal
###################################################################################
# leta <- function(DATA,t){

# }
