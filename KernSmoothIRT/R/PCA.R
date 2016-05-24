PCA<-function(x,...){

	nevalpoints   <- x$nevalpoints        # number of ticks on the x-axis
	nitem <- x$nitem	# numberof items
	evalpoints  <- x$evalpoints	# ticks on the x-axis

	X <- matrix(0,nitem,nevalpoints)       # (nitem x nevalpoints)-matrix referred to the nitem ICCs
	itemMax <- aggregate(x$OCC[,2], by = list(x$OCC[,1]),max,na.rm=TRUE)[,2]

	for(i in 1:nitem){
		

		lines0<-x$OCC[which(x$OCC[,1]==i),]


		lines1<-apply(lines0[,-c(1:3)],2,function(x)x*lines0[,3])

		
		if(x$scale[i]==0){
			
			maxitem<-max(lines0[,3])
			X[i,]<-apply(lines1,2,sum)/maxitem
		}
		else{

			X[i,]<-apply(lines1,2,sum)
		
		}
		
	}

MINX <- matrix(0,nitem,nevalpoints)
MAXX <- matrix(rep(itemMax,nevalpoints),nitem,nevalpoints)
X2   <- (X-MINX)/(MAXX-MINX)                       # normalization
R  <- apply(X2,2,rank)                                # rank 
R2 <- R - matrix(rep(colMeans(R),nitem),nitem,nevalpoints) # centering
prcomp(R2,scale=FALSE,center=FALSE,...) -> pcX
return(pcX)

}