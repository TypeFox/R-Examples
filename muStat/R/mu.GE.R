`mu.GE` <-
function(x, y=x) 
{
	if(is.array(x) && (length(dim(x))==3)) 
	{
		if((dim(x)[3]==2) && is.null(y))
			{ y <- x[,,1];  x <- x[,,2]	}
		else 
			stop("format of x and y not supported.")
	}
	else{
		if (!is.matrix(x)) x <- as.matrix(x)
		if(is.null(y)) y <- x
			else if(!is.matrix(y)) y <- as.matrix(y)
	}
    if (length(y)>1)
        apply(rbind(x,y),2,mu.GE,nrow(x))
    else 
        as.numeric(NAtoZer(outer(x[1:y],x[-(1:y)],">=")))
}
