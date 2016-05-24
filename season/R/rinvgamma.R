# rinvgamma.R from MCMCpack
rinvgamma<-function (n, shape, scale = 1) 
{
    return(1/rgamma(n, shape, scale))
}