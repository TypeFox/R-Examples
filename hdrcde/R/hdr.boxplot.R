hdr.boxplot <- function(x, prob=c(99,50), h=hdrbw(BoxCox(x,lambda),mean(prob)), lambda=1, boxlabels="", col= gray((9:1)/10), 
    main = "", xlab="",ylab="", pch=1,  ...)
{
    if(!is.list(x))
        x <- list(x)
    prob <- -sort(-prob)
    nplots <- length(x)
    junk <- unlist(x)
    ends.list <- list()
    for(j in 1:nplots)
        ends.list[[j]] <- hdr.box(x[[j]],prob,h=h,lambda=lambda)
    plot(c(0.5,nplots+0.5),range(junk,unlist(ends.list)),type="n",xlab=xlab,main=main,
        ylab=ylab,xaxt="n",...)
    if(length(boxlabels)==1 & boxlabels[1]=="")
    {
        junk <- names(x)
        if(is.null(junk))
            boxlabels <- rep("",nplots)
        else
            boxlabels <- junk
    }
    axis(1,at=c(1:nplots),labels=boxlabels,tick=FALSE,...)

    for(j in 1:nplots)
    {
        xx <- x[[j]]
        ends <- ends.list[[j]]
        nint <- length(prob)
        cols <- rep(col,5)
        for(i in 1:nint)
        {
            endsi <- ends[[i]]
            for(k in 1:(length(endsi)/2))
            {
                polygon( c(j-0.25,j-0.25,j+0.25,j+0.25,j-0.25),
                    c(endsi[k*2-1],endsi[k*2],
                    endsi[k*2],endsi[k*2-1],endsi[k*2-1]),
                    col =cols[i])
            }
        }
        for(k in 1:length(ends$mode))
            lines(c(j-0.35,j+0.35),rep(ends$mode[k],2),lty=1)
        outliers <- xx[xx<min(ends[[1]]) | xx>max(ends[[1]])]
        points(rep(j,length(outliers)),outliers,pch=pch)
    }
    invisible()
}


"hdr.box" <- function(x, prob=c(99,50), h, lambda, ...)
{
    # Does all the calculations for an HDR boxplot of x and returns
    # the endpoints of the HDR sub-intervals and the mode in a list.
    # Called by hdr.boxplot().

    r <- diff(range(x))
    if(r>0)
    {
        den <- tdensity(x,bw=h,lambda)
        info <- calc.falpha(x,den,1-prob/100)
    }
    hdrlist <- list()
    for(i in 1:length(prob))
    {
        if(r>0) hdrlist[[i]] <- hdr.ends(den,info$falpha[i])$hdr
        else hdrlist[[i]] <- c(x[1],x[1])
    }
    if(r>0) hdrlist$mode <- info$mode
    else hdrlist$mode <- x[1]
    return(hdrlist)
}

# Transformed density estimate

tdensity <- function(x,bw="SJ",lambda=1)
{
    if(is.list(x))
        x <- x[[1]]
    if(lambda==1)
        return(density(x,bw=bw,n=1001))
    else if(lambda < 0 | lambda > 1)
        stop("lambda must be in [0,1]")
    # Proceed with a Box-Cox transformed density
    y <- BoxCox(x,lambda)
    g <- density(y,bw=bw, n=1001)
    j <- g$x > 0.1 - 1/lambda # Stay away from the edge
    g$y <- g$y[j]
    g$x <- g$x[j]
    xgrid <- InvBoxCox(g$x,lambda) # x
    g$y <- c(0,g$y * xgrid^(lambda-1))
    g$x <- c(0,xgrid)
    return(g)
}

BoxCox <- function (x, lambda) 
{
    if(is.list(x))
        x <- x[[1]]
    if (lambda == 0) 
        log(x)
    else (x^lambda - 1)/lambda
}
InvBoxCox <- function (x, lambda) 
{
    if(is.list(x))
        x <- x[[1]]
    if (lambda == 0) 
        exp(x)
    else (x * lambda + 1)^(1/lambda)
}
