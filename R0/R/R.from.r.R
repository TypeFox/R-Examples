# Name   : R.from.r
# Desc   : Calculates the discretized Laplace Transform using a discretized 
#          Generation Time distribution.
# Date   : 2011/11/09
# Author : Boelle, Obadia
###############################################################################


# Function declaration

R.from.r = function#Compute the discretized Laplace Transform using a discretized Generation Time distribution
### Computes the discretized Laplace Transform using a discretized Generation Time distribution.
##details<< For internal use. Called by est.R0.EG.
##keyword<< internal

(r, ##<< exponential growth ratio
GT ##<< discretized generation time distribution
) 
  
  
# Code
  
{
	Tmax = length(GT$GT)
	R = r/sum(GT$GT * (exp(-r*(0:(Tmax-1))) - exp(-r*(1:Tmax))))
  ### An R value corresponding to inverse of discretized Laplace transform.
  ##note<< The formula for the discretized Laplace transform is taken from Wallinga.
}
