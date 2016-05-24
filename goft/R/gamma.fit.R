gammadist.fit <-  function(x)
{
    if( min(x) < 0 )
	   stop("There are negative observations. 
          \nAll data must be positive real numbers.")
    n <- length(x)
    b.check <- cov(x,log(x))
    a.check <- mean(x)/b.check
    fit <- as.matrix(c(a.check,b.check))
    colnames(fit) <- c("Parameter estimates")
    rownames(fit) <- c("shape","scale")
    return(fit)
}
