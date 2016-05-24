OLE.fun <-
function(dd, alpha){
	dd<-subset(dd, dd[,2]>0)
	sights <- rev(sort(dd[,1]))
	k <- length(sights)
	v <- (1/(k-1)) * sum(log((sights[1] - sights[k])/(sights[1] - sights[2:(k-1)])))
	e <- matrix(rep(1,k), ncol=1)
	SL<-(-log(1-alpha/2)/length(sights))^-v
	SU<-(-log(alpha/2)/length(sights))^-v
	myfun <- function(i,j,v){(gamma(2*v+i)*gamma(v+j))/(gamma(v+i)*gamma(j))}
	lambda <- outer(1:k, 1:k, myfun, v=v)
	lambda <- ifelse(lower.tri(lambda), lambda, t(lambda)) 
	a <- as.vector(solve(t(e)%*%solve(lambda)%*%e)) * solve(lambda)%*%e
	lowerCI<-max(sights) + ((max(sights)-min(sights))/(SL-1))
	upperCI<-max(sights) + ((max(sights)-min(sights))/(SU-1))
	extest<-sum(t(a)%*%sights)
	res<-data.frame(Estimate=extest, lowerCI=lowerCI, upperCI=upperCI)
	class(res)<-"extmod"
	return(res)	
	
}
