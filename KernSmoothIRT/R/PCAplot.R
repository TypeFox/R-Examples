PCAplot<-function(x,...){

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
prcomp(R2,scale=FALSE,center=FALSE) -> pcX
	#Z <- X-matrix(rep(colMeans(X),nitem),nitem,nevalpoints,byrow=T)  # Centered items (TestGraf95, p. 64)
	#Zrank <- apply(Z,2,rank) # Rank items for rank correlation
	#prcomp(Zrank,scale = F) -> pcX # Principal Component Analysis of X
	

################################
## PLOT TO ADD IN THE PACKAGE ##
################################
evalpoint <- min(evalpoints)+(max(evalpoints)-min(evalpoints))/10
par(fig=c(0.18,0.82,0.18,0.82))   
par(mai=c(0.05,0.05,0.05,0.05)) 
par(omi=c(0,0,0,0))
par(cex.axis=0.5)
plot(pcX$x[,1],pcX$x[,2],type="n",xlab="",ylab="", axes=FALSE)
box()
text(pcX$x[,1],pcX$x[,2],labels = 1:nitem,cex=0.7)
abline(v=0,lty=2)
abline(h=0,lty=2)
par(fig=c(0,0.16,0.42,0.58),new=TRUE)
plot(evalpoints,X[which(pcX$x[,1]==min(pcX$x[,1])),],type="l",xlab="",ylab="",ylim=c(0,1),axes=F)
text(evalpoint,y=0.9,labels=paste(which(pcX$x[,1]==min(pcX$x[,1])),sep=""),cex=0.6)
box()
par(fig=c(0.84,1,0.42,0.58),new=TRUE)
plot(evalpoints,X[which(pcX$x[,1]==max(pcX$x[,1])),],type="l",xlab="",ylab="",ylim=c(0,1),axes=F)
text(evalpoint,y=0.9,labels=paste(which(pcX$x[,1]==max(pcX$x[,1])),sep=""),cex=0.6)
box()
par(fig=c(0.42,0.58,0,0.16),new=TRUE)
plot(evalpoints,X[which(pcX$x[,2]==min(pcX$x[,2])),],type="l",xlab="",ylab="",ylim=c(0,1),axes=F)
text(evalpoint,y=0.9,labels=paste(which(pcX$x[,2]==min(pcX$x[,2])),sep=""),cex=0.6)
box()
par(fig=c(0.42,0.58,0.84,1),new=TRUE)
plot(evalpoints,X[which(pcX$x[,2]==max(pcX$x[,2])),],type="l",xlab="",ylab="",ylim=c(0,1),axes=F)
text(evalpoint,y=0.9,labels=paste(which(pcX$x[,2]==max(pcX$x[,2])),sep=""),cex=0.6)
box()
par(mfrow=c(1,1))
}
