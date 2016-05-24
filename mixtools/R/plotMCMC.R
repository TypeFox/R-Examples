plot.mixMCMC <- function(x, trace.plots = TRUE, summary.plots = FALSE, burnin = 2000, ...){

mix.object <- x

if (!inherits(mix.object, "mixMCMC")) 
    stop("Use only with \"mixMCMC\" objects!")

if(trace.plots==TRUE){
k<-mix.object$components
theta<-mix.object$theta
p.k=ncol(theta)
p=p.k/k
name.theta<-colnames(theta)
par(mfrow=c(p,k))
for(i in 1:p){
    for(j in 1:k){
    plot(theta[,(i-1)*k+j],type="l",ylab=name.theta[(i-1)*k+j])
    }
}
}

#regmixMH
if(is.matrix(mix.object$x) == TRUE && is.null(mix.object$y) == FALSE && summary.plots == TRUE){
y<-mix.object$y
n<-length(y)
x<-mix.object$x
p<-ncol(x)
k<-mix.object$components
theta<-mix.object$theta
if(p!=2 || sum(x[,1])!=n){                                                                  
stop(paste("This only works for simple linear regression!","\n"))
}
par(mfrow=c(1,1))
plot(x[,2],y,main="Credible Regions",xlab="Predictor",ylab="Response")
#plot(theta[-c(1:burnin),seq(1,2*k-1,by=2)],theta[-c(1:burnin),seq(2,2*k,by=2)],col=0)
for(i in 1:k){
#points(theta[-c(1:burnin),2*i-1],theta[-c(1:burnin),2*i],col=(i+1))
regcr(beta=cbind(theta[-c(1:burnin),2*i-1],theta[-c(1:burnin),2*i]),col=(i+1), x=x[,2], plot=TRUE,...)
}

}

}
