ci.RubinSchenker <-
function(x, n, conf.level) {
	z <- qnorm(1-((1-conf.level)/2))
    
    	ll <- plogis(qlogis((x+0.5)/(n+1))-z/(sqrt((n+1)*((x+0.5)/(n+1))*(1-((x+0.5)/(n+1))))))
    	ul <- plogis(qlogis((x+0.5)/(n+1))+z/(sqrt((n+1)*((x+0.5)/(n+1))*(1-((x+0.5)/(n+1)))))) 
  
    	res <- list (ci = matrix(c(ll,ul), ncol = 2))          
}
