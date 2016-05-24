opgam.iscluster.default<-function(data, idx, idxorder, alpha, ...)
{
	#Indexes of regions in the balls ordered by distance to the centre
	localidx<-idxorder[idx[idxorder]]

	localO<-sum(data$Observed[localidx])
	localE<-sum(data$Expected[localidx])
	if(localE ==0) return(c(localO, NA, NA, NA)) 

	

	localP<-sum(data$Population[localidx])

	pvalue<-ppois(localO, localE, lower.tail=FALSE)

	return (c(localO, alpha>pvalue, pvalue, sum(idx)))
}

opgam.iscluster.negbin<-function(data, idx, idxorder, alpha, mle, R=999, ...)
{
	#Indexes of regions in the balls ordered by distance to the centre
        localidx<-idxorder[idx[idxorder]]
	

	localO<-sum(data$Observed[localidx])
	localE<-sum(data$Expected[localidx])
	if(localE ==0) return(c(localO, NA, NA, NA)) 

#	bt<-boot(data[localidx, ], statistic=function(x){sum(x$Observed)},
#		R=R, sim="parametric",ran.gen=negbin.sim,
#		mle=list(n=sum(idx),size=mle$size,prob=mle$prob[localidx]))
#	pvalue<-sum(bt$t>localO)/(R+1)

	pvalue<-.Call("Ropgam_iscluster_negbin",data$Observed[localidx], data$Expected[localidx], mle$size, mle$prob[localidx], R, PACKAGE="DCluster")

	return (c(localO, alpha>pvalue, pvalue, sum(idx)))
}

opgam<-function(data, thegrid=NULL, radius=Inf, step=NULL, alpha, iscluster=opgam.iscluster.default, set.idxorder=TRUE, ...)
{
	#If the Grid is null, then create a new grid
	if(is.null(thegrid))
	{
		if(is.null(step))
			step<-.2*radius

		xgrid<-seq(min(data$x), max(data$x), by=step)
		ygrid<-seq(min(data$y), max(data$y), by=step) 

		xlen<-length(xgrid)
		ylen<-length(ygrid)
		npoints<-xlen*ylen

		thegrid<-matrix(rep(NA, 2*npoints) , ncol=2)
		
		thegrid[,1]<-rep(xgrid, times=ylen)
		thegrid[,2]<-rep(ygrid, each=xlen)
	}


	rr<-radius*radius

	GAM<-apply(thegrid, 1, opgam.intern, data, rr, set.idxorder, iscluster, alpha, ...)

	#Take only those balls which were significant
	GAM<-GAM[,!is.na(GAM[4,])]
	GAM<-as.data.frame(t(GAM[, as.logical(GAM[4,])==TRUE ]))
	#Just the first five names are set because it is possible
	#to get more than five columns if user-defined functions
	#are used
	names(GAM)<-c("x", "y", "statistic", "cluster", "pvalue", "size")

	return(GAM)
}


opgam.intern<-function(point, data, rr, set.idxorder, iscluster, alpha, ...)
{
	xd<-(data$x-point[1])
	yd<-(data$y-point[2])
	dist<-xd*xd+yd*yd

	idx<-(dist<=rr)
	if(set.idxorder) idxorder<-order(dist)

	cl<-iscluster(data=data, idx=idx, idxorder=idxorder, alpha=alpha, ...)

	return(c(point, cl))
}
