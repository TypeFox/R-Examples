getCoord <- function(p,x0,y0,nn,nb,nt){
	# Get the coordinates for canonical point (nb,nt) with mixing ratio p
	n0<-log2(nn)
	coord<-NULL
	coord$y<-log2(2*nn*(1-p)+nt*p)-1-n0+y0
	coord$x<-abs((nb*p+nn*(1-p))/(2*nn*(1-p)+nt*p)-0.5)+x0
	return(coord)
}