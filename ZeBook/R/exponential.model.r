################################################################################
# "Working with dynamic models for agriculture"
# Daniel Wallach (INRA), David Makowski (INRA), James W. Jones (U.of Florida),
# Francois Brun (ACTA)
# version : 2013-05-23
################################################################################
#' @title The Exponential growth model of dynamic of population
#' @param a : growth rate
#' @param Y0 : initial condition
#' @param duration : duration of simulation
#' @param dt : time step for integration
#' @return data.frame with Y for each time step
#' @seealso \code{\link{verhulst.update}} for the update function of the Verhulst model.
#' @export
exponential.model<-function(a,Y0,duration=40,dt=1)
{
  # Initialize variables
 # K : number of time step
 K = duration/dt+1
# 1 state variable, as 1 vectors initialized to NA
 # Y : Size of the population
 Y = rep(NA,K)
 # Initialize state variable when starting simulation on year 0
 Y[1] = Y0
# time variable
 time = 0
 # Integration loop
k=1
while (time < duration){
 	# Calculate rates of change of state variables
  	dY = a*Y[k]*dt
  	# Update state variables
  	Y[k+1]=Y[k]+dY
	# Update time
	time = time + dt
  	k=k+1
 }
  # End loop
    return(round(as.data.frame(cbind(time =seq(0,duration,by=dt) , Y)),10))
}
################################################################################
#' @title The Exponential growth model of dynamic of population - with improved Euler integration
#' @param a : growth rate
#' @param Y0 : initial condition
#' @param duration : duration of simulation
#' @param dt : time step for integration
#' @return data.frame with Y for each time step
#' @seealso \code{\link{verhulst.update}} for the update function of the Verhulst model.
#' @export
exponential.model.ie<-function(a,Y0,duration=40,dt=1)
{
# Initialize variables
# K : number of time step
 K = duration/dt+1
# 1 state variable, as 1 vectors initialized to NA
# Y : Size of the population
Y = rep(NA,K)
# Initialize state variable when starting simulation on year 0
Y[1] = Y0
# time variable
 time = 0
# Integration loop
k=1 
while (time < duration){
	# Calculate rates of change of state variables at time t
dY = a*Y[k]*dt
	# Predict value of state variables at time t+?t
YP=Y[k]+dY
# Estimate rate of change of state variables at time t+?t
dYP = a*YP*dt
    	# Predict value of state variables at time t+?t
    	Y[k+1]=Y[k]+(1/2)*(dY + dYP)
	# Update time
	time = time + dt
  	k=k+1 
  }
    # End loop
  return(round(as.data.frame(cbind(time =seq(0,duration,by=dt), Y)),10))
}
################################################################################
#' @title The Exponential growth model of dynamic of population - another form
#' @param a : growth rate
#' @param Y0 : initial condition
#' @param duration : duration of simulation
#' @param dt : time step for integration
#' @return data.frame with Y for each time step
#' @seealso \code{\link{verhulst.update}} for the update function of the Verhulst model.
#' @export
exponential.model.bis<-function(a,Y0,duration=40,dt=1)
{
	# Initialize variables
	# 1 state variable, as 1 vectors initialized to NA
	# Y : Size of the population
	Y = rep(NA,duration/dt+1)
		# Initialize state variable when starting simulation on year 0
	Y[1] = Y0
		# Integration loop
    for (k in 1:(duration/dt)){
    # Calculate rates of change of state variables  
    dY = a*Y[k]*dt
    # Update state variables 
    Y[k+1]=Y[k]+dY
  }
    # End loop
  return(round(as.data.frame(cbind(time =seq(0,duration,by=dt), Y)),10))
}
################################################################################
# end of file
