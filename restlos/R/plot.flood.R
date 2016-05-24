plot.flood <-
function (x,..., level=0)
{

# require(rgl)


if(level==0){index_drin<-x$l}else{index_drin<-which(x$fafh[,1]>=level)[1]}


#plot of input data (for d=2 and d=3) with superimposed SOM net
if(dim(x$som.results$data)[2]==2){
	dev.new()
	plot(x$som.results$data,col="red",xlab=expression(data[x]),ylab=expression(data[y]))
	points(x$som.results$code,cex=.6)
	i<-0
	repeat{i<-i+1
	lines(x$som.results$code[((i*x$som.results$xdim)-(x$som.results$xdim-1)):(i*x$som.results$xdim),],cex=.6)
	if(i==x$som.results$ydim) break}
	j<-0
	repeat{j<-j+1
	lines(x$som.results$code[c(j,(j+x$som.results$xdim)),],cex=.6)
	if(j==(x$som.results$xdim*x$som.results$ydim-x$som.results$xdim)) break}
	points(x$som.results$code[x$fafh.lib[[index_drin]],],col="GREEN",cex=.6)
	points(x$som.results$data[x$fafh.drin[[index_drin]],],col="ORANGE")
	}

if(dim(x$som.results$data)[2]==3){
	open3d(mouseMode="trackball")
	plot3d(x$som.results$data,col="red",xlab="data_x",ylab="data_y",zlab="data_z")
	points3d(x$som.results$code,size=.6)
	surface3d(x=matrix(x$som.results$code[,1],ncol=x$som.results$ydim),y=matrix(x$som.results$code[,2],ncol=x$som.results$ydim),z=matrix(x$som.results$code[,3],ncol=x$som.results$ydim),front="line",back="line")
	points3d(x$som.results$code[x$fafh.lib[[index_drin]],],col="GREEN",size=4)
	points3d(x$som.results$data[x$fafh.drin[[index_drin]],],col="ORANGE",size=4)
	}



#Mahalanobis distance plot (for data with d>3)
if(dim(x$som.results$data)[2]>3){
	dev.new()
	plot(mahalanobis(x$som.results$data,center=colMeans(x$som.results$data[x$fafh.drin[[index_drin]],]),cov=cov(x$som.results$data[x$fafh.drin[[index_drin]],])),xlab="Index",ylab="Mahalanobis distance")
	points(x$fafh.drin[[index_drin]],mahalanobis(x$som.results$data,center=colMeans(x$som.results$data[x$fafh.drin[[index_drin]],]),cov=cov(x$som.results$data[x$fafh.drin[[index_drin]],]))[x$fafh.drin[[index_drin]]],col="red")
	}



#U-landscape plot
open3d(mouseMode="trackball")
plot3d(x$umatrix,xlab="net_x",ylab="net_y",zlab="U-height")
surface3d(x=c(1:x$som.results$xdim),y=c(1:x$som.results$ydim),z=matrix(x$umatrix[,3],ncol=x$som.results$ydim),front="line",back="line")
points3d(x$umatrix[x$fafh.lib[[index_drin]],],col="RED",add=TRUE,size=4)
surface3d(x=c(1:x$som.results$xdim),y=c(1:x$som.results$ydim),z=rep(max(x$umatrix[x$fafh.lib[[index_drin]],3]),x$som.results$xdim*x$som.results$ydim),col="BLUE",add=TRUE)



#flood area flood height plot
dev.new()
plot(x$fafh[,c(1,3)],cex=.6,xlab="flood height",ylab="flood area")
lines(x$fafh[,c(1,3)],cex=.6)
if(level==0){abline(v=max(x$umatrix[x$fafh.lib[[index_drin]],3]),lty=2,lwd=2,col="red")}else{abline(v=level,lty=2,lwd=2,col="red")}


}

