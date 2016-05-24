image.hcov <-
function(x,...){	
# Initialize some variables	
Z <- x$Z
V <- x$V
diag(Z) <- 0 	
diag(V) <- 0
Z[which(abs(Z)<1e-4)] <- 0
Z[which(Z<0)] <- -1
Z[which(Z>0)] <- 1
V[which(abs(V)<1e-4)] <- 0
V[which(V<0)] <- -1
V[which(V>0)] <- 1

rgb.num=col2rgb("orange")
colororange=rgb(rgb.num[1],rgb.num[2],rgb.num[3],150,maxColorValue=256)
rgb.num=col2rgb("blue")
colorblue=rgb(rgb.num[1],rgb.num[2],rgb.num[3],100,maxColorValue=256)
colZ <- colV <- c(colorblue,"white",colororange)

if(sum(Z<0)==0){
	colZ <- c("white",colororange)
}
if(sum(V<0)==0){
	colV <- c("white",colororange)
}
if(sum(Z>0)==0){
	colZ <- c(colorblue,"white")
}
if(sum(V>0)==0){
	colV <- c(colorblue,"white")
}
if(sum(Z>0)==0 && sum(Z<0)==0){
	colZ <- c("white")
}
if(sum(V>0)==0 && sum(V<0)==0){
	colV <- c("white")
}


p <- nrow(Z)
dev.new(width=10,height=5)
set.panel(1,2)
par(oma=c(0,0,0,3))
image(t(Z),col=colZ,axes=TRUE,xaxt='n',yaxt='n',main="Z",cex.main=2,...)
image(t(V),col=colV,axes=TRUE,xaxt='n',yaxt='n',main="V",cex.main=2,...)
par(oma=c(0,0,0,1))
temp<-matrix(c(-max(abs(x$Z),abs(x$V)),max(abs(x$Z),abs(x$V)),0,0),2,2)
image.plot(temp,col=c(rep(colorblue,15),"white",c(rep(colororange,15))),horizontal=FALSE,legend.only=TRUE)		
}
