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
#       values -> array with the data (angular velocity or acceleration)
#       freq -> sampling frecuency in Herz
#
# Output Parameters:
#       data frame structure with tree fields
#         $av -> Allan Variance
#         $time -> cluster time of the computation
#         $error -> error of the estimation (quality of the variance)
#
# Description: avarn function computes the Allan variance estimation describe 
# in the equation in the documentation attached to this package.
# Therefore the input has to be the output rate/acceleration from the sensors
# In this version of the Allan variance computation the number and size of
# cluster n has been computed as the maximum number of cluster into N
# values recorded, which is ceil [(N-1)/2].
#
# License: GPL-2
#
#' @export
#
################################################################################
avarn <- function (values, freq)
{

  N = length(values)   # Number of data availables
	tau = 1/freq  # sampling time
	n = ceiling((N-1)/2) #Number of clusters

	#Allan variance array
	av <- rep (0,n)
	#Time array
	time <- rep(0,n)
	#Percentage error array of the AV estimation
	error <- rep(0,n)#1...n

  print ("Calculating...")
	# Minimal cluster size is 1 and max is n
	# in time would be 1*tau and max time would be n*tau
	for (i in 1:(n))
	{
		omega = rep(0,floor(N/i))
		T = i*tau

		l <- 1
		k <- 1
		while (k <= floor(N/i))
		{
			omega[k] = sum (values[l:(l+(i-1))])/i
			l <- l + i
			k <- k + 1
		}
		sumvalue <- 0
	
		for (k in 1: (length(omega)-1))
		{
			sumvalue = sumvalue + (omega[k+1]-omega[k])^2
		}
	
		av[i] = sumvalue/(2*(length(omega)-1))
		time[i] = T
    
		#Equation for error AV estimation
		#See Papoulis (1991) for further information
		error[i] = 1/sqrt(2*((N/i)-1))

	}

	return (data.frame(time, av, error))
}
