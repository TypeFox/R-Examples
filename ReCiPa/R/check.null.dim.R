check.null.dim <-
function(vec0){
	if(is.null(dim(vec0))){
		lov <- length(vec0)
		vec0 <- matrix(vec0,1,lov)
	}
return(vec0)
}
