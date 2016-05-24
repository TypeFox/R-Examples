loglik_joint <-
function(alpha,beta,theta,delta,x,y,R,S,family,exposure=rep(1,length(y)),negative=TRUE,zt=TRUE){
    p<-ncol(R)
    q<-ncol(S)
    n<-nrow(R)
    lambda<-as.vector(exp(S%*%beta))*exposure
	mu<-as.vector(exp(R%*%alpha))   
	#if (any(lambda<10^(-10)) | any(mu<10^(-10))) {ll<-(-10^(100))}
	#else{
	dummy<-log(density_joint(x,y,mu,delta,lambda,theta,family,zt))
    #if ((any(is.na(dummy)))| (any(dummy==-Inf))) {ll<-(-10^(100))}
	#else{
    		ll<-sum(dummy)
	#}
	#}
    if (negative==TRUE) ll<-(-ll)
    return(ll)
}
