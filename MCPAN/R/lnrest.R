lnrest <-
function(x, f)
{
    f <- as.factor(f)
    nx <- tapply(x, f, length)
    lx <- log(x)
    mlx <- tapply(lx,f,mean)
    varlx <- tapply(lx,f,var)*((nx-1)/nx)
    lmx <- mlx + 0.5 * varlx
    varlmx <- varlx/nx + (varlx^2)/(2*nx)
return(list(estimate=lmx, varest=varlmx, n=nx))}

