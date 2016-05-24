rsdepth <- function(pt, q, ...)
{
	if(!is.matrix(pt)) stop ("first argument pt must be a 2 dimensional matrix");
# check if q is a matrix of size 1x2
#	if(!is.numeric()) stop ("m must be numeric!");

	rsdepth = .C("rs_depth", 
                   xpoints=as.double(pt[,1]), 
		   ypoints=as.double(pt[,2]), 
                   query=as.double(q),  
                   depth=double(c(1)),
                   sz=as.integer(length(pt)/2),
                   PACKAGE = "rsdepth")$depth;
##print(rsdepth)

##TODO: report error in case depth is less than 0	
##	if(depth < 0) return("There was an error while calculating ");
	return(rsdepth);
}

