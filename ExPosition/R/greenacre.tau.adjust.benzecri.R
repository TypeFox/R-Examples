greenacre.tau.adjust.benzecri <-
function(eigvals,num_variables,new_eigvals,num_columns){
	Ibar <- i.bar.compute(num_variables,eigvals,num_columns)
	return((new_eigvals/Ibar)*100)
}
