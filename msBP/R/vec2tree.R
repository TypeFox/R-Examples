vec2tree <-
function(vec){
	scale <- ceiling(log2(floor(length(vec)/2)+1))
	T <- list()
	T[[1]] <- vec[1]
	j = 1
	for(s in 1:scale)
	{
	T[[s+1]] <- vec[j+1:2^s]
	j <- j + 2^s
	}
	structure(list(T=T, max.s=scale), class = "binaryTree")	
}
