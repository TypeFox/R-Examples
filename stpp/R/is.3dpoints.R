is.3dpoints <- function (x) 
{
    is <- FALSE
    	if (is.array(x)) 
      	  if (length(dim(x)) == 2) 
            	if (dim(x)[2] >= 3) 
                	is <- TRUE
    is
}
