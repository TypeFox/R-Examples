################################################################################
# "Working with dynamic models for agriculture"
# R script for practical work
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2012-04-23
# Model described in the book, Appendix. Models used as illustrative examples: description and R code
################################################################################
#' @title The Verhulst (logistic) model - calculate change for one day
#' @param Y : state variable Y(t=day)
#' @param a : growth rate
#' @param k : capacity
#' @return state variable at Y(t=day+1)
#' @seealso \code{\link{verhulst.model}} for the integration loop function of the Verhulst model.
#' @export
verhulst.update<-function(Y,a,k)
{
	# Calculate the rates of state variables (dZ)
	dY = a*Y*(1-Y/k)
  # Update the state variables Z
	Y1=Y+dY
  # Return Y1=Y(t+1)
	return(Y1)
}

################################################################################
#' @title The Verhulst (logistic) model - calculate daily values over designated time period
#' @param a : growth rate
#' @param k : capacity
#' @param Y0 : initial condition
#' @param duration : duration of simulation
#' @return data.frame with daily Y
#' @seealso \code{\link{verhulst.update}} for the update function of the Verhulst model.
#' @examples plot(verhulst.model(0.08,100,1,100), type="l", ylim=c(0,115),
#'    xlab="day", ylab="Y, population density",lwd=2)
#' @export
verhulst.model<-function(a,k,Y0,duration)
{
	# Initialize variables
	# 1 state variable, as 1 vectors initialized to NA
	# Y : Size of the population
	Y = rep(NA,duration)
		# Initialize state variable when starting simulation on year 0
	Y[1] = Y0
		# Integration loop
	for (day in 1:(duration-1))
		{
		# using the update function.
		Y[day+1] = verhulst.update(Y[day],a,k)
		}
		# End loop
  return(data.frame(day=1:duration,Y=Y[1:duration]))
}
################################################################################
# End of file
