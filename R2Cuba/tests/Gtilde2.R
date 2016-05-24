library("R2Cuba")
Gtilde2 <- function(x)
{

	cat("x", x[1], "\n")
	ret <- (x[1])^(j-1)/factorial(j-1)
        return(ret)
        
}
tt <- 3
j <- 3
NDIM <- 3
NCOMP <- 1
cuhre(NDIM, NCOMP, Gtilde2) # marche



