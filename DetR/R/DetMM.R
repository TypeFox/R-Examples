DetMM<-function(x,y,intercept=1,alpha=0.75,h=NULL,scale_est="scaleTau2",tuning.chi=1.54764,tuning.psi=4.685061){
	conv<-1
	if(!is.numeric(tuning.chi))	stop("tuning.chi nust be a finite numeric or NULL.")
	if(!is.finite(tuning.chi))	stop("tuning.chi nust be a finite numeric or NULL.")	
	if(tuning.chi<=0)		stop("tuning.chi nust be a finite numeric or NULL.")
	if(!is.numeric(tuning.psi))	stop("tuning.psi nust be a finite numeric or NULL.")
	if(!is.finite(tuning.psi))	stop("tuning.psi nust be a finite numeric or NULL.")	
	if(tuning.psi<=0)		stop("tuning.psi nust be a finite numeric or NULL.")
	s1<-robustbase::lmrob.control()
	s1$subsampling<-"simple"
	s1$max.it<-1000
	s1$k.max<-1000
	s1$maxit.scale<-1000
	niter.k<-2
	S2<-DetLTS_raw(x=x,y=y,h=h,alpha=alpha,scale_est=scale_est,doCsteps=0)
	S4<-S1<-vector("list",length(S2$Subset))
	names(S4)<-names(S1)<-names(S2$Subset)
	for(i in 1:length(which(S2$SubsetSize<nrow(x)))){
		S4[[i]]<-fast.s(x=x,y=y,int=intercept,H1=S2$Raw,k=niter.k,b=1-S2$SubsetSize[i]/nrow(x),cc=tuning.chi,conv=conv)
		names(S4[[i]])<-c('coef','scale')
		s1$bb<-1-S2$SubsetSize[i]/nrow(x)
		S1[[i]]<-robustbase::lmrob(y~x,control=s1,init=list(coefficients=S4[[i]][[1]],scale=S4[[i]][[2]]))
	}
	if(max(S2$SubsetSize)==nrow(x)){
		i<-i+1	
		S4[[i]]<-S1[[i]]<-ordreg(x=x,y=y,intercept=intercept)
	}
	S3<-vector('list',2)
	names(S3)<-c('DetS','DetMM')
	S3[[1]]<-S4
	S3[[2]]<-S1
	return(S3)
}
