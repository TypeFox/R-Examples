# Name   : plot.GT
# Desc   : A tweaked "plot" function designed to easily plot GT objects from
#          GT function
# Date   : 2011/11/09
# Author : Boelle, Obadia
###############################################################################


# Function declaration

plot.R0.GT <- function#Print the characteristics of the generation time distribution
### Prints the characteristics of the generation time distribution
##details<< For internal use. Called by print.
##keyword<< internal

(x, ##<< the generation time distribution.
... ##<< extra parameters passed to plot.
) 

  
# Code
  
{
  #Only check is that GT is of class "R0.GT"
	if (class(x)!="R0.GT") {
    stop("GT must be of class 'R0.GT'")
	}
  
  #We plot GT=f(time)
	plot(x$time,x$GT,xlab="Time",ylab="PDF",t='l', main="Generation Time distribution",...)
  
### Called for side effect. Shows a plot of the generation time distribution.	
}
