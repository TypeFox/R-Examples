################################################################################
# Automation and Robotic Section
# European Space Research and Technology Center (ESTEC)
# European Space Agency (ESA)
##
# Centre for Automation and Robotics (CAR)
# CSIC- Universidad Politecnica de Madrid
#
# Header: avarn(param1, param2)
#
# Author: Javier Hidalgo Carrio
#
# Date: 10-05-2010
#
# Input Parameters:
#       values -> array with the data recorded (angle or velocity)
#       freq -> sampling frecuency
#
# Output Parameters:
#       data frame structure with tree fields
#         $av -> Allan Variance
#         $time -> cluster time of the computation
#         $error -> error of the estimation (quality of the variance)
#
# Description: avari function computes the Allan variance estimation describe 
# in the equation in the documentation attached to this package
# Therefore the input has to be the integral value of output rate/acceleration
# from the sensors, that means angle from gyros and velocity from accelerometers
# In this version of the Allan variance computation the number and size of
# cluster n has been computed as in the avar function.
# n=2^l, l=1,2,3,...(Allan 1987)
#
# License: GPL-2
#
#' @export
#
################################################################################
avari <- function (values, freq)
{

  N = length(values) # Number of data availables
	tau = 1/freq  # sampling time
	n = ceiling((N-1)/2)
	p = floor (log10(n)/log10(2))   #Number of clusters

	#Allan variance array
	av <- rep (0,p)#0...p-1
	#Time array
	time <- rep(0,p)#0...p-1
	#Percentage error array of the AV estimation
	error <- rep(0,p)#0...p-1

  print ("Calculating...")
	# Minimal cluster size is 1 and max is 2^p
	# in time would be 1*tau and max time would be (2^p)*tau
	for (i in 0:(p-1))
	{
		omega = rep(0,N-2*(2^i)) #N-2*(2^i) is the number of that cluster
		T = (2^i)*tau
		
		sumvalue <- 0
	
		for (k in 1: (length(omega)))
		{
			sumvalue = sumvalue + (values[k+2*(2^i)]-(2*values[k+(2^i)])+values[k])^2
		}
	
		av[i+1] = sumvalue/(2*((T)^2)*(N-(2*(2^i)))) #i+1 because i starts at 0 (2^0 = 1)
		time[i+1] = T #i+1 because i starts at 0 (2^0 = 1)
		#Equation for error AV estimation
		#See Papoulis (1991) for further information
		error[i+1] = 1/sqrt(2*((N/(2^i))-1))
	}

	return (data.frame(time, av, error))
}
