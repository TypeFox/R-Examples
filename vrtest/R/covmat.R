covmat <-
function(kvec)
{
    d <- length(kvec)
    mat <- matrix(0,nrow=d,ncol=d)
    dvec <- (2*(2*kvec-1) * (kvec-1)) / (3*kvec)
    diag(mat) <- dvec

    for (i in 1:d)
    {
        for (j in 1:d)
        {
        if (i==j)
        tem <- 0
        if (j > i)
        tem <- 0
    
        mat[i,j] <- (2*(3*kvec[i]-kvec[j]-1)*(kvec[j]-1))/(3*kvec[i])
        mat[j,i] <- mat[i,j]
        }
    }
return(mat)
}
