plot.mpmcorrelogram <-
function(x,
pch=c(15,22), xlim=NULL, ylim=NULL,
        ylab=NULL, xlab=NULL,alfa=0.05,...)
        {

nclas <- length(x$clases)
  tipos <- x$pval.Bonferroni < alfa
tipos <- sapply(tipos, function(x) x=ifelse(x==TRUE,pch[1],pch[2]))
if(is.null(xlim)) xlim <- c(0,nclas +1)
if (is.null(ylim)) ylim <- c(-1,1)
if (is.null(ylab)) ylab <- expression(r[M])
if (is.null(xlab)) xlab <- "distance classes"
        plot((1:nclas)[(1:nclas)>=xlim[1] & (1:nclas)<=xlim[2]],
            x$rM[(1:nclas)>=xlim[1] & (1:nclas)<=xlim[2]],
            pch=tipos[(1:nclas)>=xlim[1] & (1:nclas)<=xlim[2]],
            cex=1.5, type="b",
            ylab=ylab, xlab=xlab, xlim= xlim,
             ylim=ylim, ...)
}

