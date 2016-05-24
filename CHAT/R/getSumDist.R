getSumDist <-
function(oo,p,sam.dat,method='both',type=1,para){
	# given a set of parameters, return the smallest summation of squared distance
	# between segment points and grid
	# the final squared distance is weighed by segment size
	# method: 'logR'--weighted by logR markers; 'BAF'--weighted by BAF markers;
	#	'both'--weighted by geometric mean of both markers.
	Ns<-dim(sam.dat)[1]
	thr<-para$thr.CL
	penalty<-para$thr.penalty
	cali<-getCanonicalLines(oo,p,type=type,para=para)
	SD<-0
	clist<-c()
	plist<-c()
	if(method=='logR'){
		w<-sam.dat[,5]
	}
	if(method=='BAF'){
		w<-sam.dat[,7]
	}
	if(method=='both'){
		w<-sqrt(sam.dat[,5]*sam.dat[,7])
	}
	nn<-type
	if(para$model==1){
		n0<-log2(1+(nn-1)*p)
		nn<-1
	}
	else{
		n0<-log2(nn)
	}
	if(para$is.tri){
		n0<-log2(2+p)-1
		b0<-abs(1/(2+p)-0.5)
		nn<-1
	}
	size.line<-100
	p0<-seq(0,1,length.out=size.line)
	nt<-1
	nb<-0
	ll<-log2(2*nn*(1-p0)+nt*p0)-1-n0+oo$y0
	bb<-abs((nb*p0+nn*(1-p0))/(2*nn*(1-p0)+nt*p0)-0.5)+oo$x0
	for(i in 1:Ns){
		if(i %in% oo$list){next}
		x<-sam.dat[i,6]
		y<-sam.dat[i,4]
		d<-c()
		for(j in 1:para$num.tracks){
			d<-c(d,Dist(x,y,cali$k[j],-1,cali$b[j]))
		}
		d.min<-which(d==min(d))
		if(d[d.min]<=thr){
			clist<-c(i,clist)
			gg<-getGrid(oo$x0,oo$y0,p,BAF=d.min-1,type=type,para=para)
			if(is.nearCP(sam.dat[i,6],sam.dat[i,4],gg,para=para)){
				plist<-c(plist,i)
			}
		}
		tmp.x<-abs(bb-x)
		delta.y<-ll[which(tmp.x==min(tmp.x))]-y
		if(delta.y>=0.5){
			# Double deletion
			d.min<-NA
		}
		if(!is.na(d.min))SD<-SD+d[d.min]^2*w[i]
	}
	SD<-SD+penalty*Dist(0,1+log(1+p)/log(2),cali$k[2],-1,cali$b[2])*(type-1)
	cali$SD<-SD/sum(w)
	cali$clist<-clist
	cali$plist<-plist
	return(cali)
}