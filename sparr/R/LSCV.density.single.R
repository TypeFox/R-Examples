LSCV.density.single <- function(h,data.ppp,res,edge,inside,ERRV=1000000){
	if(h<=0) return(NA)
	temp.dens.pts <- density(data.ppp,sigma=h,edge=edge,dimyx=c(res,res),at="points",leaveoneout=T)/data.ppp$n
	temp.dens <- density(data.ppp,sigma=h,edge=edge,dimyx=c(res,res))
	ca <- (diff(data.ppp$window$xrange)/res)*(diff(data.ppp$window$xrange)/res)
	
	if(any(is.na(temp.dens$v[t(inside)]))) return(NA)
	if(any(temp.dens$v[t(inside)]<=0)) return(ERRV)
	

	return(sum(((temp.dens$v/data.ppp$n)^2)*ca,na.rm=T)-(2/data.ppp$n)*sum(temp.dens.pts))
}