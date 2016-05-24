ci.AgrestiCoull <-
function(measure, n, conf.level) {
	z <- qnorm(1-((1-conf.level)/2))
   
   	ll <- (measure+(z^2/(2*n))-z*sqrt((measure*(1-measure)+(z^2/(4*n)))/n))/(1+(z^2/n))
   	ul <- (measure+(z^2/(2*n))+z*sqrt((measure*(1-measure)+(z^2/(4*n)))/n))/(1+(z^2/n)) 
   
   	res <- list (ci = matrix(c(ll,ul), ncol = 2))         
}
