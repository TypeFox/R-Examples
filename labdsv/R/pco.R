pco <- function(dis, k=2)
{
    tmp <-cmdscale(dis,k=k,eig=TRUE)
    class(tmp) <- "pco"
    return(tmp)
}

plot.pco <- function(x, ax = 1, ay = 2, col = 1, title = "", pch = 1, ...)
{
    if(missing(x)) {
        stop("You must specify an object of class pco from pco()")
    }
    plot(x$points[, ax], x$points[, ay], asp = 1,
        col = col, xlab = paste("PCO", ax), ylab = paste("PCO", ay),
        pch = pch, main = title, ...)
    invisible()
}


points.pco <- function(x, which, ax = 1, ay = 2, col = 2,  pch = 1, cex=1, ...)
{
    if(missing(x)) {
        stop("You must specify an object of class pco from pco()") 
    }
    if(missing(which)) {
        stop("You must specify a logical subscript")
    }
    points(x$points[, ax][which], x$points[, ay][which],col=col,pch=pch,cex=cex, ...) 
}

plotid.pco <- function(ord, ids=seq(1:nrow(ord$points)), ax = 1, ay = 2, col = 1, ...)
{
    if(missing(ord)) {
        stop("You must specify a list object from princomp()")
    }
    identify(ord$points[, ax],ord$points[, ay],ids)
}

surf.pco <- function(ord, var, ax=1, ay=2, thinplate=TRUE, col=2, labcex=0.8, 
             family=gaussian, grid=50, gamma=1.0, ...) 
{
    if (class(ord) != 'pco')
        stop("You must supply an object of class pco from pco()")
    if (missing(var)) {
        stop("You must specify a variable to surface")
    }
    x <- ord$points[,ax]
    y <- ord$points[,ay]
    if (any(is.na(var))) {
        cat("Omitting plots with missing values \n")
        x <- x[!is.na(var)]
        y <- y[!is.na(var)]
        var <- var[!is.na(var)]
    }
    if (is.logical(var)) {
        if (thinplate) {
            tmp <- gam(var~s(x,y), gamma=gamma, family=binomial)
        } else {
            tmp <- gam(var~s(x)+s(y),family=binomial,gamma=gamma)
        }
    } else {
        if (thinplate) {
            tmp <- gam(var~s(x,y), gamma=gamma, family=family)
        } else {
            tmp <- gam(var~s(x)+s(y),family=family,gamma=gamma)
        }
    }

    new.x <- seq(min(x),max(x),len=grid)
    new.y <- seq(min(y),max(y),len=grid)
    xy.hull <- chull(x,y)
    xy.hull <- c(xy.hull,xy.hull[1])
    new.xy <- expand.grid(x=new.x,y=new.y)
    inside <- as.logical(pip(new.xy$x,new.xy$y,x[xy.hull],y[xy.hull]))
    fit <- predict(tmp, type="response", newdata=as.data.frame(new.xy))
    fit[!inside] <- NA
    contour(x=new.x,y=new.y,z=matrix(fit,nrow=grid),
        add=TRUE,col=col)
    print(tmp)
    d2  <- (tmp$null.deviance-tmp$deviance)/tmp$null.deviance
    cat(paste("D^2 = ",formatC(d2,width=4),"\n"))
    invisible(tmp)
}

hilight.pco <- function (ord, overlay, ax=1, ay=2, title="", cols=c(2,3,4,5,6,7), glyph=c(1,3,5), ...)
{
    if (class(ord) != 'pco')
       stop("You must pass an object of class pco")
    overlay <- as.integer(clustify(overlay))
    plot(ord,ax=ax,ay=ay,type='n')
    title(title)
    layer <- 0
    pass <- 1
    for (i in 1:max(overlay,na.rm=TRUE)) {
        layer <- layer + 1
        if (layer > length(cols)) {
          layer <- 1
          pass <- pass + 1
        }
        col <- cols[layer]
        pch <- glyph[pass]
        points(ord, overlay == i, ax, ay, col = col, pch = pch)
    }
}

chullord.pco <- function (ord, overlay, ax = 1, ay = 2, cols=c(2,3,4,5,6,7), ltys=c(1,2,3), ...)
{
    if (class(ord) != 'pco')
        stop("You must pass an object of class pco")
    overlay <- as.integer(clustify(overlay))
    pass <- 1
    layer <- 0
    lty <- ltys[pass]
    for (i in 1:max(overlay,na.rm=TRUE)) {
        x <- ord$points[,ax][overlay==i & !is.na(overlay)]
        y <- ord$points[,ay][overlay==i & !is.na(overlay)]
        pts <- chull(x,y)
        layer <- layer + 1
        if (layer > length(cols)) {
          layer <- 1
          pass <- min(pass + 1,length(ltys))
        }
        col <- cols[layer]
        lty = ltys[pass]
        polygon(x[pts],y[pts],col=col,density=0,lty=lty)
    }
}

ellip.pco <- function (ord, overlay, ax = 1, ay = 2, 
    cols = c(2, 3, 4, 5, 6, 7), ltys = c(1, 2, 3), ...) 
{
    if (class(ord) != "pco") 
        stop("You must pass an object of class pco")
    if (inherits(overlay, c("partana", "pam", "clustering"))) {
        overlay <- overlay$clustering
        if (min(overlay) < 0 || (length(table(overlay)) != max(overlay))) {
            cat("WARNING: renumbering clusters to consecutive integers\n")
            overlay <- match(overlay, sort(unique(overlay)))
        }
    }
    else if (is.logical(overlay)) 
        overlay <- as.numeric(overlay)
    else if (is.factor(overlay)) 
        overlay <- as.numeric(overlay)
    pass <- 1
    layer <- 0
    lty <- ltys[pass]
    for (i in 1:max(overlay, na.rm = TRUE)) {
        x <- ord$points[, ax][overlay == i & !is.na(overlay)]
        y <- ord$points[, ay][overlay == i & !is.na(overlay)]
        pts <- chull(x, y)
        layer <- layer + 1
        if (layer > length(cols)) {
            layer <- 1
            pass <- min(pass + 1, length(ltys))
        }
        col <- cols[layer]
        lty = ltys[pass]
        x <- as.matrix(cbind(x[pts], y[pts]))
        elp <- ellipsoidhull(x)
        print(elp)
        lines(predict(elp),col=col)
    }
}

density.pco <- function (ord, overlay, ax = 1, ay = 2, 
    cols = c(2, 3, 4, 5, 6, 7), ltys = c(1, 2, 3), numitr=100, ...) 
{
    if (class(ord) != "pco") 
        stop("You must pass an object of class pco")
    overlay <- as.integer(clustify(overlay))
    
    densi <- function(xpts,ypts,overlay) {
        x <- xpts[overlay==1 & !is.na(overlay)]
        y <- ypts[overlay==1 & !is.na(overlay)]
        pts <- chull(x,y)
        a <- x[pts]
        a <- c(a,a[1])
        b <- y[pts]
        b <- c(b,b[1])
        inside <- pip(xpts,ypts,a,b)
        test <- pmax(inside,overlay==1)
        out <- sum(overlay)/sum(test)
        return(out)
    }

    out <- list()
    for (i in 1:max(overlay, na.rm = TRUE)) {
        obs <- densi(ord$points[, ax],ord$points[, ay], overlay==i)
        pval <- 0
        if (numitr > 0) {
            for (j in 1:(numitr-1)) {
                rnd <- sample(1:length(overlay),sum(overlay==i),replace=FALSE)
                rndvec <- rep(0,length(overlay))
                rndvec[rnd] <- 1
                tmp <- densi(ord$points[, ax],ord$points[, ay], rndvec)
                if (tmp >= obs) pval <- pval + 1
            }
        }
        pval <- (pval+1)/numitr
        print(paste('d = ',obs,'p = ',pval))
    }

}

