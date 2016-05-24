# Description: Simple Gibbs sampler demonstration on conditional exponentials from Chapter 1 (pages 30-33).
# Usage: expo.gibbs(B,k,m)
# Arguments: 	B	an upper bound
#		k	length of the subchains
#		m	number of iterations		

expo.gibbs <- function(B=5, k=15, m=5000)  {
    x <- y <- NULL
    while (length(x) < m)  {
        x.val <- c(runif(1,0,B),rep((B+1),length=k))
        y.val <- c(runif(1,0,B),rep((B+1),length=k))
        for (j in 2:(k+1))  {
            while(x.val[j] > B) x.val[j] <- rexp(1,y.val[j-1])
            while(y.val[j] > B) y.val[j] <- rexp(1,x.val[j])
        }
        x <- c(x,x.val[(k+1)])
        y <- c(y,y.val[(k+1)])
    }
    return(cbind(x,y))
}

