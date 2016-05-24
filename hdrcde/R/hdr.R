hdr <- function(x=NULL, prob=c(50,95,99), den=NULL, h=hdrbw(BoxCox(x,lambda),mean(prob)), lambda=1, nn=5000, all.modes=FALSE)
{
    if(!is.null(x))
    {
        r <- diff(range(x))
        if(r==0)
            stop("Insufficient data")
    }
    if(is.null(den))
        den <- tdensity(x,bw=h,lambda=lambda)
    alpha <- sort(1-prob/100)
    falpha <- calc.falpha(x,den,alpha,nn=nn)
    hdr.store <- matrix(NA,length(alpha),100)
    for(i in 1:length(alpha))
    {
        junk <- hdr.ends(den,falpha$falpha[i])$hdr
        if(length(junk) > 100)
        {
            junk <- junk[1:100]
            warning("Too many sub-intervals. Only the first 50 returned.")
        }
        hdr.store[i,] <- c(junk,rep(NA,100-length(junk)))
    }
    cj <- colSums(is.na(hdr.store))
    hdr.store <- matrix(hdr.store[,cj < nrow(hdr.store)],nrow=length(prob))
    rownames(hdr.store) <- paste(100*(1-alpha),"%",sep="")
    if(all.modes)
    {
        y <- c(0,den$y,0)
        n <- length(y)
        idx <- ((y[2:(n-1)] > y[1:(n-2)]) & (y[2:(n-1)] > y[3:n])) | (den$y == max(den$y))
        mode <- den$x[idx]
    }
    else
        mode <- falpha$mode
    return(list(hdr=hdr.store,mode=mode,falpha=falpha$falpha))
}

"calc.falpha" <-
function(x=NULL, den, alpha, nn=5000)
{
    # Calculates falpha needed to compute HDR of density den.
    # Also finds approximate mode.
    # Input: den = density on grid.
    #          x = independent observations on den
    #      alpha = level of HDR
    # Called by hdr.box and hdr.conf

    if(is.null(x))
        calc.falpha(x=sample(den$x, nn, replace=TRUE, prob=den$y), den, alpha)
    else
    {
        fx <- approx(den$x,den$y,xout=x,rule=2)$y
        falpha <- quantile(sort(fx), alpha)
        mode <- den$x[den$y==max(den$y)]
        return(list(falpha=falpha,mode=mode,fx=fx))
    }
}

"hdr.ends" <-
function(den,falpha)
{
    miss <- is.na(den$x) # | is.na(den$y)
    den$x <- den$x[!miss]
    den$y <- den$y[!miss]
    n <- length(den$x)
    if(falpha > max(den$y))
        return(list(falpha=falpha,hdr=NA) )
    dd <- den$y - falpha
    dd <- dd[2:n]*dd[1:(n-1)]
    index <- (1:(n-1))[dd<=0]
    index <- index[!is.na(index)]
    ni <- length(index)
    intercept <- numeric(ni)
    if(ni>0)
    {
        for(j in 1:ni)
        {
            idx <- c(index[j],index[j]+1)
            intercept[j] <- approx(den$y[idx],den$x[idx],xout=falpha)$y
        }
    }
    intercept <- sort(unique(intercept))
    ni <- length(intercept)
    if(ni == 0)
        intercept <- c(den$x[1],den$x[n])
    x1 <- 0.5*(intercept[1] + den$x[1])
    x2 <- 0.5*(intercept[ni] + den$x[n])
    if(approx(den$x,den$y,xout=x1)$y > falpha)
        intercept <- c(NA,intercept)
    if(approx(den$x,den$y,xout=x2)$y > falpha)
        intercept <- c(intercept,NA)
    return(list(falpha=falpha,hdr=intercept))
}
