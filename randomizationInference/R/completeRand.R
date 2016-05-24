#random assignment vectors according to complete randomization

completeRand=function(w,nrand){
	#formatting
	if(sum(factor(w)==w)>0) w=data.matrix(w)
	if(ncol(w)==1) w=as.vector(w)
	#random assignments
	if(is.vector(w)) lapply(1:nrand,function(i) sample(w))
	else lapply(1:nrand,function(i) w[sample(length(w[,1])),])
}