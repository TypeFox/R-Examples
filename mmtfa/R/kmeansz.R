kmeansz <-
function(x,n,G){
		kclus <- kmeans(x,G,nstart=50)$cluster
		zmat <- matrix(0,n,G)
		for(i in 1:n){
			zmat[i, kclus[i]]<-1
		}
	zmat
}
