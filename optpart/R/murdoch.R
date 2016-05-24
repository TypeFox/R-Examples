murdoch <- function (taxa,type,minval=0,minplt=10)
{
        if (!is.logical(type)) stop('The second argument must be of type logical')
        tmp <- apply(taxa>minval,2,sum)
        subtaxa <- taxa[,tmp>=minplt]
        tmp <- tmp[tmp>=minplt]
        pres <- apply(subtaxa[type,]>0,2,sum)
        abs <- apply(subtaxa[!type,]>0,2,sum)
        intype <- sum(type)
        outtype <- nrow(taxa) - intype
        first <- pres/abs
        second <- (nrow(taxa)-intype)/intype
        val <- log(first * second)
        prob <- rep(1,ncol(subtaxa))
        prob[val>=0] <- 1 - phyper(pres-1,intype,outtype,tmp)[val>=0]
        prob[val<0] <- 1 - phyper(abs-1,outtype,intype,tmp)[val<0]
        result <- list()
        result$minplt <- minplt
        result$nplots <- tmp
        result$type <- type
        result$pres <- pres
        result$abs <- abs
        result$murdoch <- val
        result$pval <- prob
        class(result) <- "murdoch"
        result
}

plot.murdoch <- function (x,axtype=1,pval=0.05,...)
{
    if (class(x) != 'murdoch') 
        stop("The first argument must be of class 'murdoch'")
    if (axtype==1) xval <- x$nplots else xval <- x$pres

    y <- x$murdoch

    if (any(y==Inf)) {
        maxy <- floor(max(y[y!=Inf])) + 2
        y[y==Inf] <- maxy
        high <- TRUE
    }
    else {
        maxy <- max(y)
        high <- FALSE
    }

    if (any(y==-Inf)) {
        miny <- ceiling(min(y[y!=-Inf])) - 2
        y[y==-Inf] <- miny
        low <- TRUE
    }
    else {
        miny <- min(y)
        low <- FALSE
    }
    labels <- pretty(range(y))
    numlbl <- as.numeric(labels)

    y[y==min(y)] <- min(numlbl)
    y[y==max(y)] <- max(numlbl)
    plot(xval,y,ylim=range(numlbl),xlab='Number of Plots',
         ylab='Murdoch Index',axes=FALSE,frame.plot=TRUE,type='n')
    if (low) {
        labels[1] <- '-Inf'
    }
    if (high) {
        labels[length(labels)] <- 'Inf'
    }
    axis(1,at=pretty(range(xval)))
    axis(2,at=pretty(range(y)),labels=labels)
    if (low) axis.break(axis=2,breakpos=(numlbl[1]+numlbl[2])/2,style='zigzag')
    if (high) axis.break(axis=2,
        breakpos=(numlbl[length(labels)]+numlbl[length(labels)-1])/2,style='zigzag')
    points(xval[y>0 & x$pval<=pval],y[y>0 & x$pval<=pval],pch="+")
    points(xval[y<0 & x$pval<=pval],y[y<0 & x$pval<=pval],pch="-")
    points(xval[x$pval>pval],y[x$pval>pval],pch=1)

    test <- readline('Do you want to identify particular species (Y or N) : ')
    if (test=='Y' || test == 'y') {
        identify(xval,y,names(x$pres))
    }
}

summary.murdoch <- function (object,pval=0.05,ndigits=3,...)
{
    if (class(object) != 'murdoch') 
        stop("The first argument must be of class 'murdoch'")
    tmp <- data.frame(round(object$murdoch,ndigits),round(object$pval,ndigits))
    tmp <- tmp[object$pval<=pval,]
    names(tmp) <- c('murdoch','pval')
    print(tmp[rev(order(tmp$murdoch)),])
}

