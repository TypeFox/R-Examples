irf <-
function(ip, x, D = 1.7){
  
	    ni = nrow(ip)
	    f = sapply(1:ni, function(i) {ip[i, 3] + (1 - ip[i, 3])/(1 + exp(-D*ip[i, 1] * (x-ip[i, 2])))})
	    r = list(x = x, f = f)
	    return(r)
	    }

