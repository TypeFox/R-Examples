mnntsplotwithmarginals<- function(cestimates,M,...)
{
	cpars<-cestimates[,3]
	R<-2
	if (length(M)!=2)
		return("Error: This is not a bivarate distribution")
	if (abs(sum(Mod(cpars)^2) - ((1/(2 * pi))^R)) > 1e-10) 
        return("Error: Sum of the squared norm of componentes greater than condition")
	x<-seq(0,2*pi-0.0000001,by=.1)
	y<-seq(0,2*pi-0.0000001,by=.1)
	m1<-max(mnntsmarginal(cestimates,M,1,x))
	m2<-max(mnntsmarginal(cestimates,M,2,y))
	m<-max(m1,m2)
	z<-outer(x,y,mnntsdensity2d,cpars,M)
	res<-persp(x,y,z,zlim=c(0,m),...)

	lines(trans3d(x,y=2*pi,z=mnntsmarginal(cestimates,M,1,x),pmat=res))
	lines(trans3d(x=0,y,z=mnntsmarginal(cestimates,M,2,y),pmat=res))

#	return(res)
}