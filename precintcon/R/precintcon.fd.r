precintcon.fd <- function(object) {
   
	if (is.element("precintcon.classified", class(object))) {
	
    ##
    # Frequency distribution
    #
		m = (object[,1] + object[,2] + 0.1)/2
		
		data <- data.frame(
			initial.class = object[,1], 
			final.class   = object[,2], 
			midpoint      = m,
			n             = object[,3],
			
			sum.n         = (function(o) {
							   
							   x <- o[1];
							   
							   for (i in 2:length(o))
								   x <- c(x, sum(o[1:i]))
								   
							   return(x)
							   
						    })(object[,3]),
						
		   P              = (function(o){
				   
							   i <- o[,1]
							   f <- o[,2]
							   n <- o[,3]
				   
							   x <- c(m[1] * n[1])								   
							   
				               for (k in 2:nrow(o))
					              x <- c(x, m[k] * n[k]) 
				   
				               return(x)
							   
			                })(object),
						
			sum.P         = (function(o){
					
							   i <- o[,1]
							   f <- o[,2]
							   n <- o[,3]
								
							   x <- c(m[1] * n[1]);								   
								
							   for (k in 2:nrow(o))
							      x <- c(x, m[k] * n[k]) 
							  
							   y <- x[1]
							   
							   for (k in 2:length(x))
							      y <- c(y, sum(x[1:k]))   
							   
							   return(y)
								
							})(object),
		     
		     p.sum.n      = (function(o) {
					 
							   x <- o[1];
								 
							   for (i in 2:length(o))
							      x <- c(x, sum(o[1:i]))
								
							   return(x/sum(o))
								 
							})(object[,3]),
						
			p.sum.P      = (function(o){
					
						      i <- o[,1]
							  f <- o[,2]
							  n <- o[,3]
								
							  x <- c(m[1] * n[1])								   
								
							  for (k in 2:nrow(o))
							     x <- c(x, m[k] * n[k]) 
							  
							  y <- x[1]/sum(x)
							 
							  for (k in 2:length(x))
							     y <- c(y, sum(x[1:k])/sum(x))
								  
							  return(y)
								
							})(object)
	    )
		
		class(data) <- c("data.frame", "precintcon.fd")
		return(data)
	
	} else {
		stop("precintcon.fd --> Object should be of class \"precintcon.classified\"");
	}
	
}


