
thull.nmds <- function (ord,var,grain,ax=1,ay=2,col=2,grid=50,nlevels=5,
                        levels=NULL,lty=1,numitr=100,...) 
{
    if (is.null(class(ord))) {
        stop("You must supply an object of class nmds from nmds")
    } else {
        if (class(ord) != "nmds") {
            stop("You must supply an object of class nmds from nmds")
        }
    }
    if(missing(var)) {
        stop("You must specify a variable to surface")
    }
    if (is.null(var)) {
        stop("No such variable")
    }
    x <- ord$points[,ax]
    y <- ord$points[,ay]
    if (any(is.na(var))) {
        cat("Omitting plots with missing values \n")
        x <- x[!is.na(var)]
        y <- y[!is.na(var)]
        var <- var[!is.na(var)]
    }
    new.x <- seq(min(x),max(x),len=grid)
    new.y <- seq(min(y),max(y),len=grid)
    hull <- matrix(0,nrow=grid,ncol=grid)

    res <- .Fortran('thull',
                    hull=as.double(hull),
                    as.double(new.x),
                    as.double(new.y),
                    as.integer(grid),
                    as.double(x),
                    as.double(y),
                    as.double(var),
                    as.integer(length(x)),
                    as.double(grain),
                    PACKAGE='labdsv')

    if (is.null(levels)) {
        vals <- levels(factor(var))
        levels <- as.numeric(vals)[-1]
    }

    contour(x=new.x,y=new.y,z=matrix(res$hull,nrow=grid),
        add=TRUE,col=col,nlevels=nlevels,levels=levels,lty=lty)
    final <- matrix(res$hull,nrow=grid)
    out <- list(thull=final,x=new.x,y=new.y,ax=x,ay=y,vals=var,
           xlab=paste('NMDS',ax),ylab=paste('NMDS',ay),
           main=deparse(substitute(var)))

    if (numitr > 0) {
        obssum <- sum(final)
        rndsum <- rep(NA,numitr-1)
        for (i in 1:(numitr-1)) {
                res <- .Fortran('thull',
                    hull=as.double(hull),
                    as.double(new.x),
                    as.double(new.y),
                    as.integer(grid),
                    as.double(x),
                    as.double(y),
                    as.double(sample(var,replace=FALSE)),
                    as.integer(length(x)),
                    as.double(grain),
                    PACKAGE='labdsv')
            rndsum[i] <- sum(res$hull) 
        } 
        cat(paste('\nvolume   = ',format(obssum,digits=5),'\nmean     = ',
                   format(mean(rndsum),digit=5),
                   '\nfraction = ',format(obssum/mean(rndsum),digits=5)))
        cat(paste('\np <= ',(sum(rndsum<=obssum)+1)/numitr),'\n')
        out$obs <- obssum
        out$reps <- rndsum
    }
    class(out) <- 'thull'
    invisible(out)
}

thull.pco <- function (ord,var,grain,ax=1,ay=2,col=2,grid=50,nlevels=5,
                       levels=NULL,lty=1,numitr=100,...) 
{
    if (is.null(class(ord))) {
        stop("You must supply an object of class pco from pco")
    } else {
        if (class(ord) != "pco") {
            stop("You must supply an object of class pco from pco")
        }
    }
    if(missing(var)) {
        stop("You must specify a variable to surface")
    }
    if (is.null(var)) {
        stop("No such variable")
    }
    x <- ord$points[,ax]
    y <- ord$points[,ay]
    if (any(is.na(var))) {
        cat("Omitting plots with missing values \n")
        x <- x[!is.na(var)]
        y <- y[!is.na(var)]
        var <- var[!is.na(var)]
    }
    new.x <- seq(min(x),max(x),len=grid)
    new.y <- seq(min(y),max(y),len=grid)
    hull <- matrix(0,nrow=grid,ncol=grid)

    res <- .Fortran('thull',
                    hull=as.double(hull),
                    as.double(new.x),
                    as.double(new.y),
                    as.integer(grid),
                    as.double(x),
                    as.double(y),
                    as.double(var),
                    as.integer(length(x)),
                    as.double(grain),
                    PACKAGE='labdsv')

    if (is.null(levels)) {
        vals <- levels(factor(var))
        levels <- as.numeric(vals)[-1]
    }

    contour(x=new.x,y=new.y,z=matrix(res$hull,nrow=grid),
        add=TRUE,col=col,nlevels=nlevels,levels=levels,lty=lty)
    final <- matrix(res$hull,nrow=grid)
    out <- list(thull=final,x=new.x,y=new.y,ax=x,ay=y,vals=var,
           xlab=paste('NMDS',ax),ylab=paste('NMDS',ay),
           main=deparse(substitute(var)))

    if (numitr > 0) {
        obssum <- sum(final)
        rndsum <- rep(NA,numitr-1)
        for (i in 1:(numitr-1)) {
                res <- .Fortran('thull',
                    hull=as.double(hull),
                    as.double(new.x),
                    as.double(new.y),
                    as.integer(grid),
                    as.double(x),
                    as.double(y),
                    as.double(sample(var,replace=FALSE)),
                    as.integer(length(x)),
                    as.double(grain),
                    PACKAGE='labdsv')
            rndsum[i] <- sum(res$hull) 
        } 
        cat(paste('\nvolume   = ',format(obssum,digits=5),'\nmean     = ',format(mean(rndsum),digit=5),
                   '\nfraction = ',format(obssum/mean(rndsum),digits=5)))
        cat(paste('\np <= ',(sum(rndsum<=obssum)+1)/numitr),'\n')
        out$obs <- obssum
        out$reps <- rndsum
    }
    class(out) <- 'thull'
    invisible(out)
}

plot.thull <- function (x,col=rainbow(20),levels=NULL,cont=TRUE,
          xlab=x$xlab,ylab=x$ylab,main=x$main,...) 
{
    if (is.null(levels)) {
        vals <- levels(factor(x$vals))
        levels <- as.numeric(vals)[-1]
    }
    image(x$x,x$y,x$thull,col=col,asp=1,xlab=xlab,ylab=ylab,main=main)
    if (cont)
        contour(x$x,x$y,x$thull,levels=levels,nlevels=length(levels),add=TRUE)
}
