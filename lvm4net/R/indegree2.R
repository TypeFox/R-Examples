indegree2 <- function(y){
	
	stopifnot(is.adjacency(y))
	
	indegree <- numeric(nrow(y))
	tabcs <- table(colSums(y))
	not0 <-  as.numeric(names(tabcs))
	indegree[not0 + 1] <- tabcs
	
	indegree
}