mh1.unif <- function( f,x0 = 0.5,x.range = c(0.01,0.99),ratio=0.1,iter = 10,... )
{   if(x.range[1] < 0.01) x.range[1] <- 0.01
    if(x.range[2] > 0.99) x.range[2] <- 0.99 
    
    chain <- rep(0,iter)
    len <- x.range[2] -x.range[1]
    chain[1] <- x0
    for(i in (1:iter)+1)
    {
	chain[i] <- chain[i-1]
	new <- runif(1,chain[i]- len * ratio,chain[i] + len * ratio)

	if(new > x.range[1] & new < x.range[2]) {
	   if(f(new,...) - f(chain[i],...) > log(runif(1)))
	      chain[i] <- new
	}	
    }
    
    chain[iter]
}

