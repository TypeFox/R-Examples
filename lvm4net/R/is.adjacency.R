is.adjacency <- function(y){
	
	is.matrix(y) && nrow(y) == ncol(y) && diag(y) == 0 && all(y %in% c(0, 1, NA))

}