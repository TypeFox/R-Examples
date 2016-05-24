
# Description: Simple Metropolis algorithm demonstration using a bivariate exponential target from Chapter 1 (pages 33-37).
# Usage: expo.metrop(m,x,y,L1,L2,L,B)
# Arguments: 	m	number of iterations		
#		x	starting point for the x vector
#		y	starting point for the y vector
# 		L1	event intensity for the x dimension
#		L2 	event intensity fot the y dimension
#		L	shared event intensity
#		B	upper bound
#		max.x	maximum for the x vector
#		max.y	maximum for the y vector

biv.exp <- function(x,y,L1,L2,L) exp( -(L1+L)*x - (L2+L)*y -L*max(x,y) )

cand.gen <- function(max.x,max.y) c(runif(1,0,max.x),runif(1,0,max.y))

expo.metrop <- function(m=5000, x=0.5, y=0.5, L1=0.5, L2=0.1, L=0.01, B=8)  {
    for (i in 1:(m-1))  {
        cand.val <- cand.gen(B,B)
        a <- biv.exp(cand.val[1],cand.val[2],L1,L2,L) / biv.exp(x[i],y[i],L1,L2,L) 
        if (a > runif(1)) { 
            x <- c(x,cand.val[1])
            y <- c(y,cand.val[2])
        }
        else  {
            x <- c(x,x[i])
            y <- c(y,y[i])
        }
    }
    return(cbind(x,y))
}

expo.metrop()




