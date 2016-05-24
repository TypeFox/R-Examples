random.shortonly.test <- function(  n=2, k=n, x.t=1, x.l=0, x.u=x.t, max.iter=1000 )
{
###
### This function generates one random short only portfolio for the given number of assets and other
### portfolio characteristics
###
### Arguments
### n = a positive integer value which is the number of assets in the portfolio
### k = a positive integer value for the number of non-zero positions
### x.t = a numeric value for the sum of the absolute value of the allocations across all assets
### x.l = a numeric value for the lower bound on the absolute value of the allocation to an asset
### x.u = a numeric value for the upper bound on the absolute value of the allocation to an asset
### max.iter = a positive integer for the maximum number of iterations in the rejection method loop
###
    thisResult <- random.longonly.test( n, k, x.t, x.l, x.u, max.iter )
    x <- -thisResult$x
    iter <- thisResult$iter
    newResult <- list( x=x, iter=iter )
    return( newResult )
}
