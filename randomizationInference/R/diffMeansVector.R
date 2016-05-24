#vector of pairwise differences of mean outcomes
#calcOptions$factors is a vector of treatment factors (defaults to a vector of 1's)
#calcOptions$pairs is a matrix of treatment level pairs (npairs x 2)

diffMeansVector=function(y,w,calcOptions){
	#formatting
	if(sum(factor(w)==w)>0) w=data.matrix(w)
	if(ncol(w)==1) w=as.vector(w)
	#calculations
	if(is.null(calcOptions)) calcOptions$pairs=matrix(c(0,1),ncol=2)
	if(is.vector(calcOptions$pairs)) calcOptions$pairs=matrix(calcOptions$pairs,ncol=2)
	if(is.null(calcOptions$factors)) calcOptions$factors=rep(1,length(calcOptions$pairs[,1]))
	sapply(1:length(calcOptions$factors),function(i)
		diffMeans(y,w,list(factor=calcOptions$factors[i],pair=calcOptions$pairs[i,])))
}
