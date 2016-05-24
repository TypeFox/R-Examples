# Name   : print.sR
# Desc   : A tweaked "print" function designed to easily print all R objects from
#          a R0.sR-class list (result of est.R0).
# Date   : 2012/04/17
# Author : Boelle, Obadia
###############################################################################

# Function declaration

print.R0.sR <- function#Plot the R0/Rt value along with confidence interval of all requested models to epidemic data
### Plot the R0/Rt value along with confidence interval of all requested models to epidemic data
##details<< For internal use. Called by print.
##keyword<< internal

(x, ##<<Result of est.R0 (class sR)) 
... ##<< Parameters passed to inner functions
 ##details<< Tweaked print() function that prints the reproduction number values for each method contained in the object constructed by est.RO().
)  
  
  # Code
  
{
  #Make sure x is of the right class.
  if (class(x)!="R0.sR") {
    stop("'x' must be of class 'R0.sR'")
  }
  
  #Successive print of individual model
  if (exists("EG", where = x$estimates)) {
    print(x$estimates$EG, ...)
  }
  
  if (exists("ML", where = x$estimates)) {
    print(x$estimates$ML, ...)
  }
  
  if (exists("AR", where = x$estimates)) {
    print(x$estimates$AR, ...)
  }
  
  if (exists("TD", where = x$estimates)) {
    print(x$estimates$TD, ...)
  }
  
  if (exists("SB", where = x$estimates)) {
    print(x$estimates$SB, ...)
  }
  
  ### Called for its side effect :
  ### Prints all R0 or R(t) values from requested estimation methods.
}
