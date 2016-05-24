subjEIS <- function(x){

	ItemExp <- matrix(NA,nrow = x$nitem, ncol=x$nsubj)

	for(it in 1:x$nitem){
		evals <-  x$OCC[which(x$OCC[,1] == it),3] %*% x$OCC[which(x$OCC[,1] == it),-c(1:3)]
		ItemExp[it,] <- approx(x=x$evalpoints, y=evals, xout=x$subjtheta)$y
	}
	
	return(ItemExp)
	
	
}