power.t <- function(n1, n2, delta, sig.level=.05){
	pt(qt(p=sig.level/2,df=n1+n2-2,lower.tail=TRUE),
	df=n1+n2-2,ncp=delta*sqrt((n1*n2)/(n1+n2)),lower.tail=TRUE)+
	pt(qt(p=sig.level/2,df=n1+n2-2,lower.tail=FALSE),df=n1+n2-2,ncp=delta*sqrt((n1*n2)/(n1+n2)),lower.tail=FALSE)
}
