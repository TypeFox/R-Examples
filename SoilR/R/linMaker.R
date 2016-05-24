#
# vim:set ff=unix expandtab ts=2 sw=2:
linMaker=function# create a linear ode system from a (solved) nonlinear one
### The function computes a (time dependent) coefficient from the pool solution and the globa solution
(
	nonLinearPoolOutput,	##<< the linearazed output operator for the pool
	globalSolution,		##<< the global solution		 
	poolSolution)		##<< the pool specific soulution
	{
	   function(Y,t){
	   Y*nonLinearPoolOutput(globalSolution(t),t)/poolSolution(t)}
### The coefficient as operator (function of Y and time). 
}
