lndest <-
function(x, f)
{
    f <- as.factor(f)
    nx <- tapply(x, f, length)
    lx <- log(x)
    mlx <- tapply(lx,f,mean)
    varlx <- tapply(lx,f,var)*((nx-1)/nx)
    mx <- exp(mlx + 0.5 * varlx)
    varmx <- (mx^2)*varlx/nx + (mx^2)*(varlx^2)/(2*nx)
return(list(estimate=mx, varest=varmx, n=nx))}

