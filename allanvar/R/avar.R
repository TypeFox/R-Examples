################################################################################
# Automation and Robotic Section
# European Space Research and Technology Center (ESTEC)
# European Space Agency (ESA)
##
# Centre for Automation and Robotics (CAR)
# CSIC- Universidad Politecnica de Madrid
#
# Header: avar(param1, param2)
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
# Description: avar function computes the Allan variance estimation describe 
# in the equation in the documentation attached to this packages.
# Therefore the input has to be the output rate/acceleration from the sensors
# In this version of the Allan variance computation the number and size of
# cluster has been computed as a power of 2 n=2^l, l=1,2,3,...(Allan 1987)
# which is convenient to estimate the amplitude of different noise components.
#
# License: GPL-2
#
#' @export
#
################################################################################
avar <- function (values, freq)
{

	N = length(values) # Number of data availables
	tau = 1/freq # sampling time
	n = ceiling((N-1)/2) 
	p = floor (log10(n)/log10(2)) #Number of clusters

	#Allan variance array
	av <- rep (0,p+1)#0 and 1...p
	#Time array
	time <- rep(0,p+1)#0 and 1...p
	#Percentage error array of the AV estimation
	error <- rep(0,p+1)#0 and 1...p

  print ("Calculating...")
	# Minimal cluster size is 1 and max is 2^p
	# in time would be 1*tau and max time would be (2^p)*tau
	for (i in 0:(p))
	{
		omega = rep(0,floor(N/(2^i))) #floor(N/(2^i)) is the number of the cluster
		T = (2^i)*tau

		l <- 1
		k <- 1
		#Perfome the average values
		while (k <= floor(N/(2^i)))
		{
			omega[k] = sum (values[l:(l+((2^i)-1))])/(2^i)
			l <- l + (2^i)
			k <- k + 1
		}
		sumvalue <- 0
	
	   #Perfome the difference of the average values
		for (k in 1: (length(omega)-1))
		{
			sumvalue = sumvalue + (omega[k+1]-omega[k])^2
		}
	
    #Compute the final step for Allan Variance estimation
		av[i+1] = sumvalue/(2*(length(omega)-1)) #i+1 because i starts at 0 (2^0 = 1)
		time[i+1] = T #i+1 because i starts at 0 (2^0 = 1)
    
		#Equation for error AV estimation
		#See Papoulis (1991) for further information
		error[i+1] = 1/sqrt(2*((N/(2^i))-1))
	}

	return (data.frame(time, av, error))
}
