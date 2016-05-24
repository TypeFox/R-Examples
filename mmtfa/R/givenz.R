givenz <-
function(n,G,known){
	zmat <- matrix(0,n,G)
	for(i in 1:G){
	  zmat[known==i, i]<-1
	}
# 	for(i in 1:n){
# 		zmat[i, known[i]]<-1
# 	}
	zmat
}
