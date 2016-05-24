drawcompletegraph <- function(x, y=NULL, startcanvas=TRUE,...)
{
	len = length(x)	
	if( is.null(y) )	{
    		if(!is.matrix(x)) stop ("first argument x must be a 2 dimensional matrix");
                   x1=as.double(x[,1])
		   y=as.double(x[,2])
		   x=x1
		   len=len/2;
	}
	##Check if there is already an active graph
	if(startcanvas)
	{
		plot(x,y)
	}
	
	for(i in 1:(length(x)) ) {
		for(j in (i+1):(length(y)) ) {
			segments(x[i],y[i],x[j],y[j],...)
		}
	}
}

