centroid <- function(x, y=NULL,...)
{
	
# check if q is a matrix of size 1x2
#	if(!is.numeric()) stop ("m must be numeric!");
	len = length(x)
	
	if(is.null(y))
	{
    		if(!is.matrix(x)) stop ("first argument pt must be a 2 dimensional matrix");
                   x1=as.double(x[,1])
		   y=as.double(x[,2])
		   x=x1
		   len=len/2;
	}

	centre = .C("rs_centroid", 
                   x=as.double(x), 
		   y=as.double(y), 
                   count=as.integer(len),  
                   cent=double(c(2)),
                   PACKAGE = c("rsdepth"))$cent;
##print(rsdepth)

##TODO: report error in case depth is less than 0	
##	if(depth < 0) return("There was an error while calculating ");
	return(centre);
}

