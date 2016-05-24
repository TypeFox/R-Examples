epolicy_loss <-
function(mu,delta,lambda,theta,family,y.max=300,zt=TRUE,compute.var=FALSE){
    if (length(mu)!= length(lambda)) stop("mu and lambda do not have the same length!")
    n<-length(mu)
    mean<-vector(length=n)
    var=NULL
    if (compute.var==TRUE){
    var=vector(length=n)
    }
    for (i in 1:n){
    #cat(paste("--- observation no ",i," ---\n"))
	foo<-function(l){
	dpolicy_loss(l,mu[i],delta,lambda[i],theta,family,y.max,zt)*l
	}
	 mean[i]<-integrate(foo,0,Inf)$value
        if (compute.var==TRUE){
	goo<-function(l){
	dpolicy_loss(l,mu[i],delta,lambda[i],theta,family,y.max,zt)*l^2
	}
	 var[i]<-integrate(foo,0,Inf)$value -mean[i]^2
        }
    }
    out=list(mean=mean,var=var)
    return(out)
}


