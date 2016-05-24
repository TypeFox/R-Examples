covmat.from.cfm <- function(Lambda, Psi.Square, tol.det=.00001)
{
# Lambda is the vector of population factor loadings
# Psi.Square is the vector of population error variances

if(length(Lambda)!=length(Psi.Square)) stop("\'Lambda\' and \'Psi.Square\' should be of the same length!")

i <- mean(length(Lambda), length(Psi.Square))

    for(r in 1:i)
    {
        if(r==1)
        {
        True.Covariance <- rep(NA, i)
        Error.Covariance <- matrix(0, i, i)
        }

    True.Covariance <- rbind(True.Covariance, Lambda[r]*Lambda)

        if(r==i)
        {
        True.Covariance <- True.Covariance[-1,]
        diag(Error.Covariance) <- Psi.Square
        Population.Covariance <- True.Covariance+Error.Covariance
        }
    }

det.value <- determinant(Population.Covariance, logarithm=FALSE)$modulus

if(det.value < tol.det) stop(paste("The determinant is less than your specified tolerance; your matrix might be singular or close to it. The determinant is ", det.value))
    
return(list(Population.Covariance=Population.Covariance, True.Covariance=True.Covariance, Error.Covariance=Error.Covariance))
}
