uniformz <-
function(n,G,clas,kno,known){
	zmat <- matrix(0, n, G)
	if(clas>0){
		for(i in 1:n){
			if(kno[i]==1){
				zmat[i, known[i]] <- 1
			}
			else{
				zmat[i,]<-1/G
			}
		}
	} 
	else{
		zmat <- zmat + 1/G
		zmat[1,] <- 0
		zmat[1,1] <- 1
	}
	zmat
}
