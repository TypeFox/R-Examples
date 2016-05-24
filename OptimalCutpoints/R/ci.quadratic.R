ci.quadratic <-
function(x, y, accuracy.measure, conf.level) {  
	if ((any (x <= 5)) | (any(y <= 5))) {
		warning(paste(accuracy.measure, " CI: \"Quadratic\" method may not be valid for some values (see Help Manual).\n", sep = ""), call. = FALSE, immediate. = TRUE)	 	
    }
	z <- qnorm(1-((1-conf.level)/2))
    
    ll <- (1/(x+y+z^2))*((x-0.5)+(z^2/2)-z*sqrt(z^2/4+((x-0.5)*(y-0.5))/(x+y)))
    ul <- (1/(x+y+z^2))*((x+0.5)+(z^2/2)+z*sqrt(z^2/4+((x+0.5)*(y+0.5))/(x+y)))
       
    res <- list (ci = matrix(c(ll,ul), ncol = 2))           
}
