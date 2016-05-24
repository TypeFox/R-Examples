#constant, additive treatment effect for factors at equally-spaced levels
#treatment indicators in 0s and 1s
#tau can be a vector

constEffect=function(y,w,w_new,poOptions){
	#formatting
	if(sum(factor(w)==w)>0) w=data.matrix(w)
	if(ncol(w)==1) w=as.vector(w)
	#calculations
	if(is.vector(w)) y+t(t((w_new-w))*poOptions$tau)
	else y+rowSums(t(t(w_new-w)*poOptions$tau))
}
