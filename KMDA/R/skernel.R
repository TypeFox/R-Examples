skernel <-
function(x,y,rho){
	if(length(x)!=length(y)) {
		cat("Error: dimensions do not match","\n")
		return (-1)
	}
	l=length(x)
	bi.x=bi.y=rep(NA,l)
	for(i in 1:l){
		bi.x[i]=as.numeric(x[i]>0)
		bi.y[i]=as.numeric(y[i]>0)
		}
	if(sum(bi.x==bi.y)<l) {out=0}   ## x and y are from different stratums
	if(sum(bi.x==bi.y)==l) out=exp(-sum(x-y)^2/rho)
	return(out)
}
