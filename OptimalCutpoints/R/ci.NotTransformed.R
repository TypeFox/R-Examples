ci.NotTransformed <-
function (x, y, measure, n, conf.level) {
	z <- qnorm(1-((1-conf.level)/2))
    
    	ll <- measure-z*sqrt((x*(1-x))/(n$d*(y^2))+(x^2*(1-y)*y)/(n$h*(y)^4))
    	ul <- measure+z*sqrt((x*(1-x))/(n$d*(y^2))+(x^2*(1-y)*y)/(n$h*(y)^4))

    	res <- list (ci = matrix(c(ll,ul), ncol = 2))    
}
