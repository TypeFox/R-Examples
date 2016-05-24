tif <-
function (ip, x, D = 1.7) 
	{
	    i = iif(ip, x, D)
	    if (is.null(dim(i$f))) 
	        dim(i$f) = c(length(i$x), length(i$f))
	    f = apply(i$f, 1, sum)
	    r = list(x = i$x, f = f, ni = ncol(i$f))
	    return(r)
	}

