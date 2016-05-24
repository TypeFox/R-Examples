# Description: 	Implementation of BCP function
# Usage:	bcp(theta.matrix,y,a,b,g,d)



bcp <- function(theta.matrix,y,a,b,g,d)  {
		    n <- length(y)
		    k.prob <- rep(0,length=n)
		    for (i in 2:nrow(theta.matrix))  {
		        lambda <- rgamma(1,a+sum(y[1:theta.matrix[(i-1),3]]), 
		                           b+theta.matrix[(i-1),3])
		        phi    <- rgamma(1,g+sum(y[theta.matrix[(i-1),3]:n]), 
		                           d+length(y)-theta.matrix[(i-1),3])
		        for (j in 1:n)  k.prob[j] <- exp(j*(phi-lambda))*
		                                     (lambda/phi)^sum(y[1:j])
		        k.prob <- k.prob/sum(k.prob)
		        k      <- sample(1:n,size=1,prob=k.prob)
		        theta.matrix[i,] <- c(lambda,phi,k)
		    }
		    return(theta.matrix)
		} 

