compRev <-
function(vec){
	return(as.array(apply(apply(matrix(comp(apply(vec,1,s2c)),nchar(vec[1]), length(vec)),2,rev),2,c2s)))
}
