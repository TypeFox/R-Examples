ci.wald <-
function(x, y, accuracy.measure, measure, n, conf.level) {  
	if ((any (x <= 20)) | (any(y <= 20))) {
      	warning(paste(accuracy.measure, " CI: \"Wald\" method may not be valid for some values (see Help Manual).\n", sep = ""), call. = FALSE, immediate. = TRUE)	
	}
    z <- qnorm(1-((1-conf.level)/2))        
    
    ll <- measure-(z*sqrt((measure*(1-measure))/n)+1/(2*n))
    ul <- measure+(z*sqrt((measure*(1-measure))/n)+1/(2*n))
    
    res <- list (ci = matrix(c(ll,ul), ncol = 2))      
}
