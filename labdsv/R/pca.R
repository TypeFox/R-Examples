pca <- function(mat, cor=FALSE, dim=min(nrow(mat),ncol(mat)))
{
    tmp <- prcomp(mat, retx=TRUE, center=TRUE, scale=cor)
    out <- list()
    out$scores <- tmp$x[,1:dim]
    out$loadings <- tmp$rotation[,1:dim]
    out$sdev <- tmp$sdev[1:dim]
    out$totdev <- sum(tmp$sdev^2)
    class(out) <- "pca"
    return(out)
}
                

plot.pca <- function(x, ax = 1, ay = 2, col = 1, title = "", pch = 1, ...)
{
    if (class(x) != 'pca')
        stop("You must specify a an object of class pca")
    plot(x$scores[, ax], x$scores[, ay], asp=1, 
        col = col, xlab = paste("PCA", ax), ylab = paste("PCA", ay), 
        pch = pch, main = title)
    invisible()
}

points.pca <- function(x, which, ax = 1, ay = 2, col = 2,  pch = 1, cex = 1, ...)
{
    if (class(x) != 'pca')
        stop("You must specify a list object from pca")
    if (missing(which)) {
        stop("You must specify a logical subscript")
    }
    if (length(which) != nrow(x$scores)) {
        stop("Points specifier must be of the same length as the number of samples")
    }
    points(x$scores[, ax][which], x$scores[, ay][which],col=col,pch=pch,cex=cex) 
}

plotid.pca <- function(ord, ids=seq(1:nrow(ord$scores)), ax = 1,  ay = 2, col = 1, ...)
{
    if (class(ord) != 'pca')
        stop("You must specify a list object from princomp()")
    identify(ord$scores[, ax],ord$scores[, ay],ids)
}

surf.pca <- function(ord, var, ax=1, ay=2, col=2, labcex=0.8, 
                   family=gaussian, thinplate=TRUE, grid=50, gamma=1.0, ...)
{
    if (class(ord) != 'pca') 
        stop("You must specify an object of class pca")
    if(missing(var)) {
        stop("You must specify a variable to surface")
    }
    x <- ord$scores[, ax]
    y <- ord$scores[, ay]
    if (is.logical(var)) {
        if (thinplate) {
            tmp <- gam(var~s(x,y),family=binomial,gamma=gamma)
        } else {
            tmp <- gam(var~s(x)+s(y),family=binomial,gamma=gamma)
        }
    } else {
        if (thinplate) {
            tmp <- gam(var~s(x,y),family=family,gamma=gamma)
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


summary.pca <- function(object, dim=length(object$sdev), ...)
{
    vars <- object$sdev^2
    vars <- vars/object$totdev
    cat("Importance of components:\n")
    print(rbind("Standard deviation" = object$sdev[1:dim],
        "Proportion of Variance" = vars[1:dim],
        "Cumulative Proportion" = cumsum(vars[1:dim])))
}

scores.pca <- function (x,labels=NULL,dim=length(x$sdev)) 
{
    if (dim>length(x$sdev)) {
        cat("Only",length(x$sdev)," axes available\n")
        dim <- length(x$sdev)
    }
    if (!is.null(labels)) {
        cbind(labels,x$scores[,1:dim])
    } else {
        x$scores[,1:dim]
    }
}

loadings.pca <- function (x, dim=length(x$sdev), digits=3, cutoff=0.1)
{
    if (dim>ncol(x$loadings)) {
        cat("Only",ncol(x$loadings),"axes available\n")
        dim <- ncol(x$loadings)
    }
    cat("\nLoadings:\n")
    cx <- format(round(x$loadings[,1:dim], digits = digits))
    cx[abs(x$loadings[,1:dim]) < cutoff] <- substring("       ",1, nchar(cx[1, 1]))
    print(cx, quote = FALSE)
    invisible()
}

varplot.pca <- function(x,dim=length(x$sdev)) 
{
    var <- x$sdev^2
    barplot(var[1:dim],ylab="Variance")
    readline("Hit Return to Continue\n")
    barplot(cumsum(var/x$totdev)[1:dim],ylab="Cumulative Variance")
}

hilight.pca <- function (ord, overlay, ax=1, ay=2, cols=c(2,3,4,5,6,7), glyph=c(1,3,5), ...)
{
    if (class(ord) != 'pca')
       stop("You must pass an object of class pca")
    overlay <- as.integer(clustify(overlay))
    plot(ord,ax-ax,ay=ay,type='n')
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

chullord.pca <- function (ord, overlay, ax = 1, ay = 2, cols=c(2,3,4,5,6,7), ltys=c(1,2,3),
...)
{
    if (class(ord) != 'pca')
        stop("You must pass an object of class pca")
    overlay <- as.integer(clustify(overlay))
    pass <- 1
    layer <- 0
    lty <- ltys[pass]
    for (i in 1:max(overlay,na.rm=TRUE)) {
        x <- ord$scores[,ax][overlay==i & !is.na(overlay)]
        y <- ord$scores[,ay][overlay==i & !is.na(overlay)]
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

