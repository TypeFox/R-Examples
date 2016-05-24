

Kappa<-function(CMX,st.dev=TRUE){
### CMX is a confusion matix from which the function will calculate Kappa (proportion of 
### specific agreement) and its standard deviation.
###
###	CMX		confusion matrix from cmx()
###	st.dev	should the standard deviation be calculated

### check data format

if(nrow(CMX)!=ncol(CMX) || is.matrix(CMX)==FALSE){
	stop("'CMX' must be a confusion matrix")}

### check logicals

if(is.logical(st.dev)==FALSE){
	stop("'st.dev' must be of logical type!")}

### Check for NA values

if(sum(is.na(CMX))!=0){return(NA)}

### Do calculations

	N<-sum(CMX)
	N.cat<-length(CMX[1,])
	R.sum<-apply(CMX,1,sum)
	C.sum<-apply(CMX,2,sum)
	Diag<-diag(CMX)
	theta.1<-sum(diag(CMX))/N
	theta.2<-sum(R.sum*C.sum)/(N^2)
	KAPPA<-(theta.1-theta.2)/(1-theta.2)
	if(st.dev==FALSE){
		return(Kappa=KAPPA)
	}else{
		theta.3<-sum(Diag*(R.sum+C.sum))/(N^2)
		theta.4<-0
		for(i in 1:N.cat){
			for(j in 1:N.cat){
				theta.4<-theta.4+(CMX[i,j]*((C.sum[i]+R.sum[j])^2))}}
		theta.4<-theta.4/(N^3)
	
		theta.1b<-(1-theta.1)
		theta.2b<-(1-theta.2)
		KV.1<-(theta.1*theta.1b)/(theta.2b^2)
		KV.2<-(2*theta.1b*((2*theta.1*theta.2)-theta.3))/(theta.2b^3)
		KV.3<-((theta.1b^2)*(theta.4-(4*(theta.2^2))))/(theta.2b^4)
		KAPPA.var<-((KV.1+KV.2+KV.3)/N)
		eps <- 10^(-16)
		if(is.finite(KAPPA.var)==TRUE){
			if(abs(KAPPA.var)<eps){KAPPA.var<-0}}
		KAPPA.sd<-KAPPA.var^.5
		return(data.frame(Kappa=KAPPA,Kappa.sd=KAPPA.sd))}
}
