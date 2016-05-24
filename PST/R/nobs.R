## number of observations

setMethod("nobs", "PSTf", function(object) {
	
	res <- sum(object[[1]][["e"]]@n)
	
	return(res)
})




