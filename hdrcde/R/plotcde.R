plot.cde <- function(x, firstvar=1, mfrow=n2mfrow(dim(x$z)[firstvar]), plot.fn="stacked",x.name,margin=NULL,...)
{
    dimz <- dim(x$z)
    if(plot.fn=="hdr")
        fn <- hdr.cde
    else if(plot.fn=="stacked")
        fn <- stacked.plot
    else
        stop("Unknown plotting function")

    if(length(dimz)<2)
        stop("Insufficient data")
    else if(length(dimz)==2)
        return(invisible(fn(x,...)))
    else if(length(dimz)>3)
        stop("More than 2 explanatory variables not yet implemented")
    else if(length(dimz)-1 != length(x$x))
        stop("Dimension of z does not match x")

    if(missing(x.name))
        x.name <- x$x.name
    oldpar <- par(mfrow=mfrow)
    small.den <- x
    index <- firstvar + (firstvar==1) - (firstvar==2)
    small.den$x <- x$x[[index]]
    small.den$x.name <- x.name[index]
    xlabels <- x$x[[firstvar]]
    if(max(xlabels)-min(xlabels) > 30)
        xlabels <- round(xlabels)
    else if(max(xlabels) - min(xlabels) > 3)
        xlabels <- round(xlabels,1)
    else if(max(xlabels) - min(xlabels) > 0.3)
        xlabels <- round(xlabels,2)
    for(i in 1:dimz[firstvar])
    {
        if(firstvar==1)
            small.den$z <- x$z[i,,]
        else
            small.den$z <- x$z[,i,]
        if(!is.null(margin))
            mden <- margin[,i]
        else
            mden <- rep(1,length(small.den$x))
        fn(small.den,main=paste(x.name[firstvar],"=",xlabels[i]),mden=mden,...)
    }
    par(mfrow=oldpar$mfrow)
    invisible()
}

stacked.plot <- function(den,mden=rep(1,length(den$x)),threshold=0.05,main="",xlab,ylab,font=1,col=4,view=c(4,1,3),...)
{
    if(missing(ylab))
        ylab <- den$y.name
    if(missing(xlab))
        xlab <- den$x.name

    # Omit missing values
    if(!is.matrix(den$z))
        den$z <- matrix(den$z,nrow=1)
    nz <- dim(den$z)
    miss <- (1:length(den$y))[is.na(apply(den$z,length(nz),sum))]
    if(length(miss)>0)
    {
        den$y <- den$y[-miss]
        if(length(nz)==2)
            den$z <- den$z[,-miss]
        else if(length(nz)==3)
            den$z <- den$z[,,-miss]
        else
            stop("I can't handle more than 2 covariates")
    }


    if(length(den$x)==1)
    {
        main <- paste(ylab,"|",xlab,"=",den$x)
        plot(den$y,den$z[1,],type='l',xlab=ylab,ylab="Density")
        mtext(main,3,cex=1,font=font)
        return(invisible())
    }

    oldpar <- par()
    par(mar=c(1,0,0,1)+.1)
    maxden <- max(den$z,na.rm=TRUE)
    zlim <- c(0,maxden)
    xlim <- range(unlist(den$x),na.rm=TRUE)
    ylim <- range(den$y,na.rm=TRUE)
    xrange <- xlim[2]-xlim[1]
    yrange <- ylim[2]-ylim[1]
    junk <- persp(xlim,ylim,cbind(rep(0,2),rep(maxden,2)), zlim=zlim,box=FALSE,axes=FALSE,
        theta=-75,phi=25,border=NA,...)
    mtext(main,3,cex=1,font=font,line=-4)


    # Function to enable things to be added to perspective plot
    perspp <- function(x,y,z, pmat)
    {
      tr <- cbind(x,y,z,1) %*% pmat
      list(x = tr[,1]/tr[,4], y= tr[,2]/tr[,4])
    }

    # add axes
    lines(perspp(xlim,c(ylim[1],ylim[1]),c(0,0),junk))
    lines(perspp(c(xlim[1],xlim[1]),ylim,c(0,0),junk))
    lines(perspp(xlim,c(ylim[2],ylim[2]),c(0,0),junk))
    lines(perspp(c(xlim[2],xlim[2]),ylim,c(0,0),junk))
    # add tick marks and scale
    ylabels <- pretty(ylim)
    ylabels <- ylabels[ylabels<=max(ylim) & ylabels >= min(ylim)]
    for(i in ylabels)
    {
        lines(perspp(c(xlim[1],xlim[1] - xrange/75),c(i,i),c(0,0),junk))
        lines(perspp(c(xlim[1],xlim[1] + xrange/75),c(i,i),c(0,0),junk))
        text(perspp(xlim[1]-xrange/20,i,0,junk),paste(round(i,4)),cex=0.8,font=font)
    }
    xlabels <- pretty(xlim)
    xlabels <- xlabels[xlabels<=max(xlim) & xlabels >= min(xlim)]
    for(i in xlabels)
    {
        lines(perspp(c(i,i),c(ylim[1],ylim[1]-yrange/75),c(0,0),junk))
        lines(perspp(c(i,i),c(ylim[2],ylim[2]+yrange/75),c(0,0),junk))
        text(perspp(i,ylim[1]-yrange/20,0,junk),paste(round(i,4)),cex=0.8,font=font)
    }
    midx <- 0.5*(xlim[2]+xlim[1])
    midy <- 0.5*(ylim[2]+ylim[1])
    text(perspp(midx,ylim[1] - yrange/7, 0,junk),xlab,cex=0.8,srt=265,font=font)
    text(perspp(xlim[1] - xrange/7,midy,0,junk),ylab, cex=0.8,srt=352,font=font)
    n <- length(den$x)
    if(view[1]>0)
        index <- n:1
    else
        index <- 1:n
    cols <- rep(col,n)[1:n]
    maxden <- max(mden)
    for(i in index)
    {
        if(length(c(den$z[i,])[!is.na(den$z[i,])]) > 0 & mden[i] > threshold*maxden)
        {
            pol <- perspp(rep(den$x[i],length(den$y)+2),
            c(den$y[1],den$y,den$y[length(den$y)]),
            c(0,.5*den$z[i,],0),junk)
            polygon(pol$x, pol$y, col=cols[i], border=TRUE)
        }
    }
    par(mar=oldpar$mar)
    par(new=FALSE)
    invisible(junk)
}

hdr.cde <- function(den,prob=c(50,95,99),plot=TRUE,plot.modes=TRUE,
            mden=rep(1,length(den$x)),threshold=0.05,nn=1000,
            xlim,ylim,xlab,ylab, border=TRUE,font=1,cex=1,...)
{
    if(missing(ylab))
        ylab <- den$y.name
    if(missing(xlab))
        xlab <- den$x.name
    if(missing(cex))
        cex=1
    nx <- nrow(den$z)
    na <- length(prob)
    hdrr <- list()
    modes <- matrix(NA,nx,4)
    for(i in 1:nx)
    {
        sumz <- sum(den$z[i,],na.rm=TRUE)
        if(!is.na(sumz))
        {
            if(sumz > 0)
            {
                info <- hdr(den=list(x=den$y,y=den$z[i,]),prob=prob,nn=nn)
                hdrr[[i]] <- t(info$hdr)
                modes[i,] <- c(info$mode,rep(NA,4-length(info$mode)))
            }
        }
    }
    if(!plot)
        return(list(hdr=hdrr,modes=modes))
    nalpha <- ncol(hdrr[[1]])
    nint <- length(hdrr)
    x.margin <- den$x
    if(missing(ylim))
        ylim <- range(unlist(hdrr),na.rm=TRUE)
    if(length(x.margin)>1)
        xwidth <- x.margin[2]-x.margin[1]
    else
        xwidth <- 1e-5
    if(missing(xlim))
        xlim <- range(x.margin)+c(-xwidth,xwidth)
    plot.default(x.margin,x.margin*0,type="n",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,font=font,cex=cex,xaxt="n",...)
    if(length(x.margin)>1)
    {
        axis(1)
        axis(3,labels=FALSE)
    }
    else
        axis(1,at=x.margin)
    axis(4, labels = FALSE)
    cols <- gray((9:1)/10)
    maxden <- max(mden)
    index <- rep(FALSE,nint)
    for(j in 1:nint)
    {
        # test if density high enough to plot HDR
        if(mden[j] > threshold*maxden)
        {
            for(k in 1:nalpha)
                add.hdr(hdrr[[j]][,k], x.margin[1] + (j-1)*xwidth, xwidth,
                        border=border, col=cols[k+1])
            index[j] <- TRUE
        }
    }
    if(plot.modes)
        matpoints(x.margin[index],matrix(modes[index,],nrow=sum(index)),cex=cex,col=1,pch=19)
    invisible(return(list(hdr=hdrr,modes=modes)))
}
