
com.ci <- function( b, varb, crit, z, n, l1, l2, method, lambda, res.df, r.factor)
{
    b.ci <- c(NA,NA)
    arguments <- list( l1=l1, l2=l2, z=z, n=n, method=method, crit=crit, lambda=lambda, res.df=res.df, r.factor=r.factor )
    #check if limits have opposite sign
    val.b <- lr.b.com(b,arguments)
    b.m   <- b - 2*sqrt(crit*varb)
    b.p   <- b + 2*sqrt(crit*varb)
    val.m <- lr.b.com(b.m,arguments)
    val.p <- lr.b.com(b.p,arguments)
    #if necessary, move limits further from b (small sample issues)
    if ( val.m*val.b > 0 )
        {
         b.m <- b - 4*sqrt(crit*varb)
         val.m <- lr.b.com(b.m,arguments)
        }
    if ( val.p*val.b > 0 )
        {
         b.p <- b + 4*sqrt(crit*varb)
         val.p <- lr.b.com(b.p,arguments)
        }
    #if still problems, move further and adjust for possibly 0 variance
    if ( val.m*val.b > 0 )
        {
         b.m <- b - 8*sqrt(crit*(varb+0.1))
         val.m <- lr.b.com(b.m,arguments)
        }
    if ( val.p*val.b > 0 )
        {
         b.p <- b + 8*sqrt(crit*(varb+0.1))
         val.p <- lr.b.com(b.p,arguments)
        }
    res <- uniroot( lr.b.com, c(b.m, b ), tol = 0.0001, arguments=arguments )
    b.ci[1] <- res$root
    res <- uniroot( lr.b.com, c(b , b.p ), tol = 0.0001, arguments=arguments )
    b.ci[2] <- res$root

    if ( b.ci[1]==b.ci[2] )
    {
        str("Same limits - unable to find a different solution!?")
    }
    b.ci
}
