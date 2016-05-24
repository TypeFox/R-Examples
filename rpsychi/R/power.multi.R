power.multi <- function(n, n.ind, delta, sig.level=.05){
	v <- n-n.ind-1
	lamda <- delta*(n.ind+v+1)
	pf(qf(p=sig.level, df1=n.ind, df2=v,lower.tail=FALSE),df1=n.ind, df2=v, ncp=lamda, lower.tail=FALSE)
}
