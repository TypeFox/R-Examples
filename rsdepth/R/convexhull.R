convexhull <- function(x, y=NULL,...)
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
	
	
	id=chull(x,y)
	sz=length(id)
	x2=matrix(rnorm(sz))
	y2=matrix(rnorm(sz))

	for(i in 1:sz)
	{

	x2[i] = x[id[i]]
	y2[i] = y[id[i]]

	}
	return(cbind(x2,y2));
}

