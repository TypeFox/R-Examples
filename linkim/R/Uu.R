Uu <-
function(A){ 
	A <- data.frame(A)
	c <- length(A)
	A <- unlist(A)
	A <- matrix(A,,c)#dim(A)
	return(A)
}
