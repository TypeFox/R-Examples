getOrigin <- function(sam.dat,ntry=10,concensus.thr=0.6,para){
	thr<-para$thr.originsize
	origin<-NULL
	origin.cluster<-rep(0,dim(sam.dat)[1])
	dat.norm <- (range(sam.dat[,4])[2]-range(sam.dat[,4])[1])/(range(sam.dat[,6])[2]-range(sam.dat[,6])[1])
	o.list<-c()
	k.thr<-para$thr.kmeans
	if(para$is.seed){
		while(1){
			for(k in 1:ntry){
				cl<-getKmeans(sam.dat,para=para)
				ncen<-length(cl$centers)/2
				centers<-cl$centers
				min.dist<-999
				zeroP <- -1
				for(i in 1:ncen){
					o.size<-sqrt(sum(sam.dat[cl$cluster==i,5]*sam.dat[cl$cluster==i,7]))
					dist <- sqrt(centers[i,1]^2+(centers[i,2]/dat.norm)^2)
					if(dist<=min.dist & o.size>=para$thr.originsize){
						min.dist<-dist
						zeroP<-i
					}
				}
				origin.cluster[cl$cluster==zeroP]<-origin.cluster[cl$cluster==zeroP]+1
			}
			o.list<-which(origin.cluster>=ntry*concensus.thr)
			if(length(o.list)>0){break}
			k.thr<-k.thr*2
		}
		x0<-weighted.mean(sam.dat[o.list,6],sam.dat[o.list,7])
		y0<-weighted.mean(sam.dat[o.list,4],sam.dat[o.list,5])
		origin.cluster[which(origin.cluster<ntry*concensus.thr)]<-0
		origin.cluster[o.list]<-1
	}
	else{
		x0<-0
		y0<-0
	}
	while(1){
		tag<-0
		for(i in 1:length(origin.cluster)){
			if(origin.cluster[i]==1){next}
			x<-sam.dat[i,6]
			y<-sam.dat[i,4]
			if(abs(x-x0)<=para$std.BAF & abs(y-y0)<=para$std.LRR){
				tag<-1
				origin.cluster[i]<-1
				o.list<-c(o.list,i)
				if(length(o.list)>1){
					x0<-weighted.mean(sam.dat[o.list,6],sam.dat[o.list,7])
					y0<-weighted.mean(sam.dat[o.list,4],sam.dat[o.list,5])
				}
				else{
					x0<-x/2
					y0<-y/2
				}
			}
		}
		if(tag==0){break}
	}
	origin$list<-o.list
	origin$x0<-x0
	origin$y0<-y0
	return(origin)
}
