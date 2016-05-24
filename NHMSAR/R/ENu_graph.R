ENu_graph <-
function(data,u,lty=1,col=1,add=FALSE,CI=FALSE,alpha=0.05,N.s.data=NULL,xlab = "P(Y<u)",ylab="Intensity of upcrossings",ylim=NULL){
F = NULL
T = dim(data)[1]
N.samples = dim(data)[2]
d <- dim(data)[3]
if(is.null(d)|is.na(d)){d <- 1}
data = array(data,c(T,N.samples,d))
Nu = matrix(0,length(u),N.samples)
CI.mat = matrix(0,length(u),2)

for (iu in 1:length(u)) {
	F[iu] = sum(data<u[iu])
	for (ex in 1:N.samples) {
		Nu[iu,ex] = sum(data[1:(T-1),ex,1]>u[iu] & data[2:T,ex,1]<u[iu] )
	}
	if (CI) {
		tmp = matrix(Nu[iu,],N.s.data,N.samples/N.s.data)
		CI.mat[iu,] = quantile(apply(tmp,2,sum),probs=c(alpha/2,1-alpha/2))
	}
}
F = F/length(data)
Nu1 = apply(Nu,1,sum)/N.samples

if (add == FALSE) {
	if (is.null(ylim)) {plot(c(F,1),c(Nu1,0),typ="l",lty=lty,col=col,xlab = xlab,ylab=ylab,lwd=3)} else
	{plot(c(F,1),c(Nu1,0),typ="l",lty=lty,col=col,xlab = xlab,ylab=ylab,lwd=3,ylim=ylim)}
}
else {lines(c(F,1),c(Nu1,0),lty=lty,col=col,lwd=3/2)}
if (CI) {
	CI.mat = CI.mat/N.s.data
	lines(c(F,1),c(CI.mat[,1],0),lty=3,col=col)
	lines(c(F,1),c(CI.mat[,2],0),lty=3,col=col)
}
list(u=u,F=F,Nu=Nu1,CI = CI.mat)
}
