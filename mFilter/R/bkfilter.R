### Baxter-King filter
bkfilter <- function(x,pl=NULL,pu=NULL,nfix=NULL,type=c("fixed","variable"),drift=FALSE)
{
    if(is.null(drift)) drift <- FALSE
    xname=deparse(substitute(x))
    type = match.arg(type)

    if(is.null(type)) type <- "fixed"

    if(is.ts(x))
        freq=frequency(x)
    else
        freq=1

    if(is.null(pl))
    {
        if(freq > 1)
            pl=trunc(freq*1.5)
        else
            pl=2
    }

    if(is.null(pu))
        pu=trunc(freq*8)

    b = 2*pi/pl
    a = 2*pi/pu

    n = length(x)

    if(n<5)
        warning("# of observations in Baxter-King filter < 5")

    if(pu<=pl)
        stop("pu must be larger than pl")
    if(pl<2)
    {
        warning("in Baxter-King kfilter, pl less than 2 , reset to 2")
        pl = 2
    }

    if(is.null(nfix))
        nfix = freq*3

    if(nfix>=n/2)
        stop("fixed lag length must be < n/2")

    j = 1:(2*n)
    B = as.matrix(c((b-a)/pi,(sin(j*b)-sin(j*a))/(j*pi)))

    AA = matrix(0,n,n)

    if(type=="fixed")
    {
        bb = matrix(0,2*nfix+1,1)
        bb[(nfix+1):(2*nfix+1)] = B[1:(nfix+1)]
        bb[nfix:1] = B[2:(nfix+1)]
        bb = bb-sum(bb)/(2*nfix+1)

        for(i in (nfix+1):(n-nfix))
            AA[i,(i-nfix):(i+nfix)] = t(bb)
    }

    if(type=="variable")
    {
        for(i in (nfix+1):(n-nfix))
        {
            j=min(c(i-1,n-i))
            bb=matrix(0,2*j+1,1)
            bb[(j+1):(2*j+1)] = B[1:(j+1)]
            bb[j:1] = B[2:(j+1)]
            bb = bb-sum(bb)/(2*j+1)
            AA[i,(i-j):(i+j)] = t(bb)
        }
    }
    xo = x
    x = as.matrix(x)
    if(drift)
        x = undrift(x)

    x.cycle = AA%*%as.matrix(x)
    x.cycle[c(1:nfix,(n-nfix+1):n)] = NA
    x.trend = x-x.cycle
    if(is.ts(xo))
    {
        tsp.x = tsp(xo)
        x.cycle=ts(x.cycle,star=tsp.x[1],frequency=tsp.x[3])
        x.trend=ts(x.trend,star=tsp.x[1],frequency=tsp.x[3])
        x=ts(x,star=tsp.x[1],frequency=tsp.x[3])
    }
    res <- list(cycle=x.cycle,trend=x.trend,fmatrix=AA,title="Baxter-King Filter",
                xname=xname,call=as.call(match.call()),
                type=type,pl=pl,pu=pu,nfix=nfix,method="bkfilter",x=x)

    return(structure(res,class="mFilter"))
}
