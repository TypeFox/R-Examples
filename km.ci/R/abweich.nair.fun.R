"abweich.nair.fun" <-
function(matrix,c)  
{
    # Calculates the upper and lower derivation to the Kaplan-Meier estimator
    # for determining the boundaries of a Hall-Wellner band (which is
    # symmetric).
    
    kap.mei <- matrix[,2]
    sigma <- matrix[,3]
    result1 <- c*sqrt(sigma)*kap.mei
    result2 <- exp((c*sqrt(sigma))/log(kap.mei))
    return(list(lin.dev=result1,log.dev=result2))
}

