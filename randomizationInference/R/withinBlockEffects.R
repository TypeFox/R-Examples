#vector of within-block pairwise differences of mean outcomes
#calcOptions$block is a vector of block designations
#calcOptions$factors is a vector of treatment factors (defaults to a vector of 1's)
#calcOptions$pairs is a matrix of treatment level pairs (npairs x 2)
#calcOptions$blockindex is a vector of block indices (defaults to a vector of 1's)

withinBlockEffects=function(y,w,calcOptions){
	#formatting
	if(sum(factor(w)==w)>0) w=data.matrix(w)
	if(ncol(w)==1) w=as.vector(w)
	#calculations
	if(is.null(calcOptions)) calcOptions$pairs=matrix(c(0,1),ncol=2)
	if(is.vector(calcOptions$pairs)) calcOptions$pairs=matrix(calcOptions$pairs,ncol=2)
	if(is.null(calcOptions$factors)) calcOptions$factors=rep(1,length(calcOptions$pairs[,1]))
	if(is.null(calcOptions$blockindex)) calcOptions$blockindex=rep(1,,length(calcOptions$pairs[,1]))
	sapply(1:length(calcOptions$factors),function(i)
		diffMeans(y[calcOptions$block==calcOptions$blockindex[i]],
		w[calcOptions$block==calcOptions$blockindex[i]],
		list(factor=calcOptions$factors[i],pair=calcOptions$pairs[i,])))
}