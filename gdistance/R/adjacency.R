# Author: Jacob van Etten jacobvanetten@yahoo.com
# Date :  November 2013
# Version 1.2
# Licence GPL v3

adjacencyFromTransition <- function(x)
{
	tc <- transitionCells(x)
	x <- transitionMatrix(x)
	transition.dgT <- as(x,"dgTMatrix")
	adj <- cbind(transition.dgT@i+1,transition.dgT@j+1)
	adj <- cbind(tc[adj[,1]], tc[adj[,2]])
	return(adj)
}