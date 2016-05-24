SubstDesign <- function(PX2 = NULL,PX1 = NULL){

	if(is.null(PX1) | is.null(PX2)  ){stop("The user must enter two design sets")}

	d <- dim(as.matrix(PX2))[2]
	n2 <- dim(as.matrix(PX2))[1]
	n1 <- dim(as.matrix(PX1))[1]
	
	dist <- 0 
	for(i in 1:d){
		grid <- expand.grid(PX2[,i],PX1[,i])
		dist <- dist + (grid[,1]-grid[,2])^2
	}
	Matdist <- matrix(dist,n2,n1)
	indice <- max.col(-(Matdist))
	
	PX1 <- as.matrix(PX1[-indice,])
	PX1 <- rbind(PX1,PX2)

	return(list(PX = PX1,le = length(indice)))
}