run_ppp <- function(data,xy,h,WIN,counts){
	ppp.full <- ppp(x=data[,1],y=data[,2],window=WIN,check=F)
	ppp.den.QEX <- density(x=ppp.full,sigma=h,xy=xy,weights=counts*rep(1/sum(counts),nrow(data)),spill=1)
	return(ppp.den.QEX)
}