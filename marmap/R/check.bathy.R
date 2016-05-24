check.bathy = function(x){
	order(as.numeric(colnames(x))) -> xc
	order(as.numeric(rownames(x))) -> xr
	x[xr, xc] -> sorted.x
	return(sorted.x)
}