ci.transformed <-
function (x, y, n, conf.level) {
	z <- qnorm(1-((1-conf.level)/2))
     
    	ll <- exp(log(x/y)-z*sqrt((1-x)/(n$d*x)+(1-y)/(n$h*y)))
    	ul <- exp(log(x/y)+z*sqrt((1-x)/(n$d*x)+(1-y)/(n$h*y))) 
    
    	res <- list (ci = matrix(c(ll,ul), ncol = 2))     
}
