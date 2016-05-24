PQln <-
function(nx, mlx, varlx, B=10000)
{
    Zx<-rnorm(n=B, mean=0, sd=1)
    Chix<-rchisq(n=B, df=nx-1)
    Tx <- mlx - (Zx*sqrt(varlx))/((sqrt(Chix)/sqrt(nx-1))*sqrt(nx)) + (varlx)/(2*Chix/(nx-1))
    return(Tx)
}

