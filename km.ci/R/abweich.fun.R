"abweich.fun" <-
function(matrix,k,n)  
{
    # Calculates the upper and lower derivation to the Kaplan-Meier estimator
    # for determining the boundaries of a Hall-Wellner band (which is
    # symmetric).
    
    kap.mei <- matrix[,2]
    sigma <- matrix[,3]
    result1 <- k*(1+n*sigma)*kap.mei/sqrt(n)
    result2 <- exp(k*(1+n*sigma)/(sqrt(n)*log(kap.mei)))
    return(list(lin.dev=result1,log.dev=result2))
}

