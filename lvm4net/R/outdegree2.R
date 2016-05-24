outdegree2 <- function(y){
	
	stopifnot(is.adjacency(y))
	
	outdegree <- numeric(nrow(y))
	tabrs <- table(rowSums(y))
	not0 <-  as.numeric(names(tabrs))
	outdegree[not0 + 1] <- tabrs
	
	outdegree
	
}