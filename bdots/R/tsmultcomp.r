tsmultcomp = function(rho, alpha, N, df, maxItr=1000, tol=1e-8,verbose=FALSE)
{
    LB <- 0
    UB <- 1
    x0 <- alpha
    x1 <- effectiveAlpha(rho, x0, N, df)
    diff <- abs(x1-alpha)
    itrs <- 0
    while((diff > tol || x1 > alpha) && itrs < maxItr)
    {   
        if (x1 > alpha)
        {
            UB <- x0
            x0 <- (LB+x0)/2
        }
        else if(x1 < alpha)
        {
            LB <- x0
            x0 <- (UB+x0)/2
        }
        if (verbose)
        {
            cat(paste("Iteration: ", itrs, "\n", sep = ""))
        }
        x1 <- effectiveAlpha(rho, x0, N, df)
        diff <- abs(x1-alpha)

        itrs <- itrs + 1
    }
		if(x0 > .05) x0 <- .05
    return(x0)
}