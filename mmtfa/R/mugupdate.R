mugupdate <-
function(G,zmat,w,x,p,mug,n){
	for(g in 1:G){
#   zw <- 
#		mug[g,] <- colSums(c(zmat[,g])*c(w[,g])*x)/sum(c(zmat[,g])*c(w[,g]))
    
#		mug[g,] <- colSums(zmat[,g]*w[,g]*x)/sum(zmat[,g]*w[,g])
#	  zw <- zmat[,g]*w[,g]
#	  mug[g,] <- .colSums(zw*as.matrix(x),n,p)/sum(zw)
	  mug[g,] <- .colSums(zmat[,g]*w[,g]*as.matrix(x),n,p)/sum(zmat[,g]*w[,g])
	}
	mug
}
