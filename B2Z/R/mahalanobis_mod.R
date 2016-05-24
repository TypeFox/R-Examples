##################################################
#This function computes the mahalanobis distance #
#with respect to a matrix 'invcov'               #
##################################################

mahalanobis_mod <- function (x, d, center, invcov) 
{
    x <- matrix(x, ncol = d)
    x <- sweep(x, 2, center)
    retval <- rowSums((x %*% invcov) * x)
    return(retval)
}