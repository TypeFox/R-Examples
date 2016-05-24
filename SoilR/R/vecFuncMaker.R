#
# vim:set ff=unix expandtab ts=2 sw=2:
vecFuncMaker=function# creates a vector valued function from the functions for the components
### The function is a helper to create a vector valued function of two arguments which is very useful to create systems of ode
(
funcs,	##<< The list of functions computing the vector components
arg1,   ##<< The first argument of the component functions
arg2	##<< The second argument of the component functions
)
{
	function(arg1,arg2){
		matrix(byrow=TRUE,
			nrow=length(funcs),
			mapply(
				function(fun){fun(arg1,arg2)},
				funcs
			)
		)
	}
	### A vector valued function with the vector size equal to the number of 
	### functions in the first argument
}
