MLE <-
function(resp, ip, D = 1.7, min = -4, max = 4)
	{
		np = nrow(resp)
		
		logf<-function (x, r, p) {
	    	pr = p[,3] + (1 - p[,3])/(1 + exp(-D*p[, 1] * (x-p[, 2])))
	    	ll = r * log(pr) + (1 - r) * log(1 - pr)
	    	lf = sum(ll)
	    	return(lf)} 
		
		esti<-function(x,resp,ip) optimize(logf, lower = min, upper = max, maximum = TRUE, r = resp, p = ip)$maximum
	
	  o = sapply(1:np, function(i) esti(resp = resp[i, ],ip=ip))
	 
	  return(o)
	 }

