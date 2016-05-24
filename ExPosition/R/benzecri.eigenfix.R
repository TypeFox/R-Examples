benzecri.eigenfix <-
function(eigvals,num_variables){
	new_eigvals <- eigvals[eigvals > (1/num_variables)]
	new_eigvals <- (num_variables/(num_variables-1)) * (new_eigvals -(1/num_variables))
	new_eigvals <- new_eigvals^2
	return(new_eigvals[new_eigvals > (2*.Machine$double.eps)])
}
