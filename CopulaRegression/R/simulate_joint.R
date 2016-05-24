simulate_joint <-
function(n,mu,delta,lambda,theta,family,max.y=5000,eps=0.00001,zt=TRUE){
    # simulate gamma-distributed observations
    x<-rgam(n,mu,delta)
    y<-vector(length=n)
    for (i in 1:n){
    stop=FALSE
    upper=max(50,floor(lambda))
    while (stop==FALSE){
    # compute probabilities for 1:upper
    if (zt==TRUE){
        yy=1:upper
    }
    if (zt==FALSE){
        yy=0:upper
    }
    my.prob<-density_conditional(x=x[i],y=yy,mu,delta=delta,lambda=lambda,theta=theta,family=family,zt=zt)
    # check if the probability for Y=upper is small enough or if upper is larger than the maximal number 
    if (upper>max.y |my.prob[upper]<=eps) {stop=TRUE} else {upper=upper*5}
    }
    # check that probabilities are not zero
    non.zero<-(my.prob!=0)
    y[i]<-sample(yy[non.zero],1,prob=my.prob[non.zero])
    }
    Z=cbind(x,y)
    colnames(Z)=c("x","y")
    return(Z)
}
