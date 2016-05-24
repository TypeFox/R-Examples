cv <- function(C.of.V=NULL, mean=NULL, sd=NULL, N=NULL, unbiased=FALSE)
{

if(!is.null(mean) & !is.null(sd))
{
k <- sd/mean
if(!is.null(C.of.V)) stop("Since \'mean\' and \'sd\' were specified, do not specify \'C.of.V\'.")
}

if(unbiased==TRUE)
{    if (!is.null(C.of.V))
    k <- C.of.V* (1+1/(4*N))
    if (is.null(C.of.V))
    k <- k*(1+1/(4*N))
}

return(k)
}
