sigma.hat<-function(y,f1,f2){
	if (length(f2)>1){
		z<-matrix(nrow=dim(y)[1],ncol=dim(y)[2])
		n1<-nlevels(as.factor(f1))
		n2<-nlevels(as.factor(f2))
		for (i in 1:n1){
			for (j in 1:n2){
				ix<-which(f1==i&f2==j)
				if (length(ix)>1){
					v<-apply(y[,ix],1,mean)
				}
				if (length(ix)==1){
					v<-y[,ix]
				}
				z[,ix]<-y[,ix]-v
			}
		}
		return(z%*%t(z)/(length(f1)-n1*n2))
	}
	if (length(f2)==1){
		z<-matrix(nrow=dim(y)[1],ncol=dim(y)[2])
		n1<-nlevels(as.factor(f1))
		for (i in 1:n1){
			ix<-which(f1==i)
			if (length(ix)>1){
				v<-apply(y[,ix],1,mean)
			}
			if (length(ix)==1){
				v<-y[,ix]
			}
			z[,ix]<-y[,ix]-v
		}
		return(z%*%t(z)/(length(f1)-n1))
}
}