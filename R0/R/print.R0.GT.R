# Name   : print.GT
# Desc   : A tweaked "print" function designed to print useful data on GT objects 
#          from GT function
# Date   : 2011/11/09
# Author : Boelle, Obadia
###############################################################################


# Function declaration

print.R0.GT <- function#Print the characteristics of the generation time distribution
### Prints the characteristics of the generation time distribution.
##details<< For internal use. Called by print.
##keyword<< internal

(x, ##<< the GT distribution.
... ##<< Parameters passed to inner functions
 ) 
  
  
# Code

{
  #Make sure x is of the right class
	if (class(x) != "R0.GT") {
    stop("x must be of class 'R0.GT'")
  }
	
  cat("Discretized Generation Time distribution\n")
	cat("mean:",x$mean,", sd:",x$sd,"\n")
	print.default(x$GT)
	cat("\n")
  
  ### Called for side effect. Displays GT and mean/sd.
}
