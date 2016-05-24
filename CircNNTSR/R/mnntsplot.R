#mnntsplot<- function(cpars,M,...)
mnntsplot<- function(cestimates,M,...)
{
	cpars<-cestimates[,3]
	R<-2
	if (length(M)!=2)
		return("Error: This is not a bivarate distribution")
	if (abs(sum(Mod(cpars)^2) - ((1/(2 * pi))^R)) > 1e-10) 
        return("Error: Sum of the squared norm of componentes greater than condition")
	x<-seq(0,2*pi-0.0000001,by=.1)
	y<-seq(0,2*pi-0.0000001,by=.1)
	z<-outer(x,y,mnntsdensity2d,cpars,M)
	res<-persp(x,y,z,...)
#	return(res)
}