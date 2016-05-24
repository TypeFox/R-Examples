### Butterworth filter
bwfilter <- function(x,freq=NULL,nfix=NULL,drift=FALSE)
{

    if(is.null(drift)) drift <- TRUE
    xname=deparse(substitute(x))
    if(is.ts(x))
        frx=frequency(x)
    else
        frx=1

    if(is.null(freq))
    {
        if(frx > 1)
            freq=trunc(frx*2.5)
        else
            freq=2
    }

    if(is.null(nfix))
        nfix = 2

    xo = x
    x = as.matrix(x)
    if(drift)
        x = undrift(x)

    n = length(x)

    cut.off = 2*pi/freq
    mu = (1/tan(cut.off/2))^(2*nfix)

    imat = diag(n)
    Ln = rbind(matrix(0,1,n),diag(1,n-1,n))
    Ln = imat-Ln
    if(nfix > 1)
    {
        for(i in 1:(nfix-1))
            Ln = (imat-Ln)%*%Ln
    }
    Q  = t(Ln[3:n,])
    SIGMA.R = t(Q)%*%Q
    SIGMA.n = abs(SIGMA.R)
    g = t(Q)%*%as.matrix(x)
    b = solve(SIGMA.n+mu*SIGMA.R,g)
    x.cycle = c(mu*Q%*%b)
        x.trend = x-x.cycle

    if(is.ts(xo))
    {
        tsp.x = tsp(xo)
        x.cycle=ts(x.cycle,star=tsp.x[1],frequency=tsp.x[3])
        x.trend=ts(x.trend,star=tsp.x[1],frequency=tsp.x[3])
        x=ts(x,star=tsp.x[1],frequency=tsp.x[3])
    }
    A = mu*Q%*%solve(SIGMA.n+mu*SIGMA.R)%*%t(Q)

        if(is.ts(xo))
    {
        tsp.x = tsp(xo)
        x.cycle=ts(x.cycle,star=tsp.x[1],frequency=tsp.x[3])
        x.trend=ts(x.trend,star=tsp.x[1],frequency=tsp.x[3])
        x=ts(x,star=tsp.x[1],frequency=tsp.x[3])
    }

    res <- list(cycle=x.cycle,trend=x.trend,fmatrix=A,title="Butterworth Filter",
                xname=xname,call=as.call(match.call()),
                type="asymmetric",lambda=mu,nfix=nfix,freq=freq,method="bwfilter",x=x)

    return(structure(res,class="mFilter"))
}
