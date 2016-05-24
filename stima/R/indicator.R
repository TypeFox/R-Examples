indicator <-
function(x){
	x<-as.factor(x)
	x<-as.double(x)
	G<-matrix(rep(0,length(x)*max(x)),ncol=max(x))
	for (i in 1:max(x)){
				G[x==i,i]<-1
		}
		return(G)
	}
