getKmeans <- function(sam.dat,iter.max=10,para){
	# BAF data are first normalized so that it has the same range as logR
	# a threshold is set to grasp the first best-fit kmeans center numbers
	# this threshold is normalized by logR range
	thr<-para$thr.kmeans
	dat.norm <- (range(sam.dat[,4])[2]-range(sam.dat[,4])[1])/(range(sam.dat[,6])[2]-range(sam.dat[,6])[1])
	thr.norm <- thr*dat.norm
	BAF<-sam.dat[,6]*dat.norm
	logR<-sam.dat[,4]
	ncenters <- 2
	while(1){
		cl <- kmeans(cbind(BAF,logR),iter.max=iter.max,centers=ncenters,nstart=50)
		tag <- 0
		for(i in cl$withinss){
			if(i >= thr.norm){
				tag<-1
				break
			}
		}
		if(tag==0){
			break
		}
		ncenters <- ncenters+1
	}
	cl$centers[,1]<-cl$centers[,1]/dat.norm
	return(cl)
}