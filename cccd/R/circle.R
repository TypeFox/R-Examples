draw.circle <- function(x,y,radius,border=NULL,col=NA,nv=100,...)
{
	angle.inc <- 2*pi/nv
   theta <- seq(0,2*pi-angle.inc,by=angle.inc)
	if(length(radius)<length(x)) radius <- rep(radius,length.out=length(x))
	if(length(col)<length(x)) col <- rep(col,length.out=length(x))
	if(!is.null(border))
		if(length(border)<length(x)) border <- rep(border,length.out=length(x))
	for(i in 1:length(x)){
	   polygon(radius[i]*cos(theta)+x[i],radius[i]*sin(theta)+y[i],
		        border=border[i],col=col[i],...)
	}
}
