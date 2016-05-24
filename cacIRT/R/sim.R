sim <-
function(ip, x, D = 1.7){
	    i = irf(ip, x, D)
    	d = dim(i$f)
   		u = runif(d[1] * d[2])
    	dim(u) = d
    	return(ifelse(i$f > u, 1, 0))}

