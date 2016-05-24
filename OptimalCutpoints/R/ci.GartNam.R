ci.GartNam <-
function(x, y, n, conf.level) {
	z <- qnorm(1-((1-conf.level)/2))
    
    	ll <- (2*n$d*n$h*x*y-sqrt((2*n$d*n$h*x*y)^2-4*(n$h*n$d*y^2-z^2*n$d*y*(1-y))*(n$h*n$d*x^2-z^2*n$h*x*(1-x))))/(2*(n$h*n$d*y^2-z^2*n$d*y*(1-y)))
    	ul <- (2*n$d*n$h*x*y+sqrt((2*n$d*n$h*x*y)^2-4*(n$h*n$d*y^2-z^2*n$d*y*(1-y))*(n$h*n$d*x^2-z^2*n$h*x*(1-x))))/(2*(n$h*n$d*y^2-z^2*n$d*y*(1-y)))
    
    	res <- list (ci = matrix(c(ll,ul), ncol = 2))    
}
