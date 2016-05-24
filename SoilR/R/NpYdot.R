#
# vim:set ff=unix expandtab ts=2 sw=2:
NpYdot=function
### All the ode solvers need to be able to compute the derivative y'=y'(y,t) for all possible values of y and t.
### This funtion is used to create such a function suitable for the use by the odesolvers available in the package.  (\code{\link{euler}} or \code{\link{ode}} or any other user provided function with the same interface.)
### The returned function y' defines the ode system and thus reflects the model.
(A	##<< An object representing a decomposition operator, basically a matrix valued function containing the decay rates for the n pools , coupling and feedback paths of the whole model as functions of time. The size of this quadratic matrix is nxn for n pools.
,
inputrates ##<< An object of Class inputFlux, basically a  vector function describing the inputs to the pools as functions of time
 ){
   ydot=function(y,t){
      D=(A(t)%*%y)+inputrates(t)
   }
   return(ydot)
   ### The return value is a function (or closure in R jargon ) that can compute the derivative vector for the ode system  y'=y'(y,t).
}
