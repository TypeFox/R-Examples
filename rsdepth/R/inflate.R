inflate <- function(x, y=NULL, factor=2, ...)
{
	len = length(x)
	
	if(is.null(y))
	{
    		if(!is.matrix(x)) stop ("first argument pt must be a 2 dimensional matrix");
                   x1=as.double(x[,1])
		   y=as.double(x[,2])
		   x=x1
		   len=len/2;
	}
	
	x2=matrix(rnorm(len))
	y2=matrix(rnorm(len))

	##plot(x*2,y*2)
	pt=centroid(x,y)
	
	for(i in 1:len)
	{
		x2[i]=((x[i]-pt[1])*factor)+pt[1];
		y2[i]=((y[i]-pt[2])*factor)+pt[2];
	}


	return(cbind(x2,y2));
}

