#
# vim:set ff=unix expandtab ts=2 sw=2:
solver=function
### This function applies one of its arguments (which is an ode solver) to its other arguments, the righthand side of an ode system and a vector of the initial values and the times where a value is sought.
### Its main purpose is to provide a very easy and flexible interface to the toplevel functions that use it.
### A user can thus decide which method is used to solve his model by providing the method as an argument.
(
 times,			##<< A vector containing the points in time where the solution is sought.
 ydot,			##<< The function used by the odesolver to compute the derivative. For an example how such a function should look like have a look at \code{\link{NpYdot}} which creates such a function for the n pool example. 
 startValues,		##<< a vector containing the initial amounts of carbon for the n pools. The length of this vector has to be equal to the number of pools and thus compatible with ydot. At the moment this is not checked by the function.

 solverFunc=deSolve.lsoda.wrapper		##<< The function used by to actually solve the ODE system. This can be \code{\link{euler}} or \code{\link{ode}} or any other user provided function with the same interface. 
 )
 {
   Y=solverFunc(times,ydot,startValues) 
   return(Y)
}
