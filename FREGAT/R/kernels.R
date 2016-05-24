# Functions from package 'SKAT' v.0.82 (c) 2011

call_Kernel_IBS<-function(Z,n,p){
	K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
	temp<-.C("Kernel_IBS", as.integer(as.vector(t(Z))), as.integer(n), as.integer(p), as.double(as.vector(K)), PACKAGE='FREGAT')[[4]]
	matrix(temp,nrow=n)
}

call_Kernel_IBS_Weight<-function(Z,n,p,weights){
	given_weight = 1
	if( is.null(weights)){
		weights = rep(0,p);	given_weight = 0
	} else {
		weights<-weights^2
	}
	K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
	temp<-.C("Kernel_IBS_Weight",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.integer(given_weight),
	as.double(weights),as.double(as.vector(K)), PACKAGE='FREGAT')[[6]]
	matrix(temp,nrow=n)
}

call_Kernel_2wayIX<-function(Z,n,p){
	K<- matrix(rep(0,n*n),nrow = n, ncol = n)	
	temp<-.C("Kernel_2wayIX",as.integer(as.vector(t(Z))),as.integer(n), as.integer(p),as.double(as.vector(K)), PACKAGE='FREGAT')[[4]]
	matrix(temp,nrow=n)
}

lskmTest.GetKernel = function(Z, kernel, weights,n,m){
    if (kernel == 'quadratic')        K = (Z%*%t(Z)+1)**2
	if (kernel == 'IBS')              K = call_Kernel_IBS(Z,n,m)
    if (kernel == 'IBS.weighted')     K = call_Kernel_IBS_Weight(Z,n,m,weights)
  	if (kernel == '2wayIX')           K = call_Kernel_2wayIX(Z,n,m)  
	return(K)
}
