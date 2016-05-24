getDistToPath <- function(x0,y0,x,y,type,nb,nt,strictness=1.1,thr=0.001,is.scaling=F,BAF.scale=8){
	# Find the nearest distance from segment (x,y) to a canonical path
	# thr: threshold for p to call for the nearest distance
	pst<-0
	ped<-strictness
	coord.st<-getCoord(pst,x0,y0,nn=type,nb,nt)
	coord.ed<-getCoord(ped,x0,y0,nn=type,nb,nt)
	if(nb==nt/2|nt==2*type)is.scaling<-FALSE
	if(is.scaling){
		#re-scale BAF so that it has the same dimension as LRR
		BAF.scale<-(coord.ed$y-coord.st$y)/(coord.ed$x-coord.st$x)
	}
	#plot(0,0,cex=0,xlim=c(coord.st$x,coord.ed$x),ylim=c(coord.st$y,coord.ed$y))
	#plot(0,0,cex=0,xlim=c(0,0.5),ylim=c(-4,4))
	#points(x,y,col=4,pch=8)
	while(ped-pst>=thr){
		p.half<-(ped+pst)/2
		coord.st<-getCoord(pst,x0,y0,nn=type,nb,nt)
		coord.ed<-getCoord(ped,x0,y0,nn=type,nb,nt)
		coord.half<-getCoord(p=p.half,x0,y0,nn=type,nb,nt)
		#points(coord.st,col=2,cex=0.1)
		#points(coord.ed,col=2,cex=0.1)
		v1<-c((coord.ed$x-coord.st$x)*BAF.scale,coord.ed$y-coord.st$y)
		v2<-c((x-coord.half$x)*BAF.scale,y-coord.half$y)
		if(getCosine(v1,v2)<0)ped<-p.half
		else{pst<-p.half}
	}
	minDist<-list(p.min=pst,dist=sum(c((coord.st$x-x)*BAF.scale,coord.st$y-y)^2),x.min=coord.st$x,y.min=coord.st$y)
	return(minDist)
}