i.bar.compute <-
function(num_variables,eigvals,num_columns){
	return((num_variables/(num_variables-1)) * ((sum(eigvals) - (num_columns-num_variables)/(num_variables^2))))
}
