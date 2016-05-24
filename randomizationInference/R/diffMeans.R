#pairwise difference of mean outcomes
#input: outcomes (y), assignments (w),
#	calcOptions=list of options for calculating test statistic (if necessary),
#output: value of test statistic
#calcOptions$factor can denote the treatment factor to evaluate (1 is default)
#calcOptions$pair can denote the pair of treatment levels to compare (c(0,1) is default)

diffMeans=function(y,w,calcOptions=NULL){
	#formatting
	if(sum(factor(w)==w)>0) w=data.matrix(w)
	if(ncol(w)==1) w=as.vector(w)
	#calculations
	if(is.null(calcOptions$pair)) calcOptions$pair=c(0,1)
	if(is.vector(w)){
		mean(y[w==calcOptions$pair[2]])-mean(y[w==calcOptions$pair[1]])
	}else{
		if(is.null(calcOptions$factor)) calcOptions$factor=1
		mean(y[w[,calcOptions$factor]==calcOptions$pair[2]])-mean(y[w[,calcOptions$factor]==calcOptions$pair[1]])
	}
}