#' plclust in colour
#' 
#' Modifiction of plclust for plotting hclust objects in *in colour*!
#' 
#' @param hclust hclust object
#' @param labels a character vector of labels of the leaves of the tree
#' @param lab.col colour for the labels; NA=default device foreground colour
#' @param hang as in \code{\link{hclust}} & \code{\link{plclust}}
#' @param xlab title for x-axis (defaults to no title)
#' @param sub subtitle (defualts to no subtitle)
#' @param ... further arguments passed to \code{\link{plot}}
#' @author Eva KF Chan
#' 
#' @examples
#' data(iris)
#' hc <- hclust( dist(iris[,1:4]) )
#' myplclust(hc, labels=iris$Species,lab.col=as.numeric(iris$Species))
#' 
#' 
myplclust <- function( hclust, labels=hclust$labels, lab.col=rep(1,length(hclust$labels)), hang=0.1, xlab="", sub="", ...){
 ## modifiction of plclust for plotting hclust objects *in colour*!
 ## Copyright Eva KF Chan 2009
 ## Arguments:
 ##    hclust:    hclust object
 ##    lab:        a character vector of labels of the leaves of the tree
 ##    lab.col:    colour for the labels; NA=default device foreground colour
 ##    hang:     as in hclust & plclust
 ## Side effect:
 ##    A display of hierarchical cluster with coloured leaf labels.
    y <- rep(hclust$height,2)
    x <- as.numeric(hclust$merge)
    y <- y[which(x<0)]
    x <- x[which(x<0)]
    x <- abs(x)
    y <- y[order(x)]
    x <- x[order(x)]
    plot( hclust, labels=FALSE, hang=hang, xlab=xlab, sub=sub, ... )
    text( x=x, y=y[hclust$order]-(max(hclust$height)*hang), labels=labels[hclust$order], col=lab.col[hclust$order], srt=90, adj=c(1,0.5), xpd=NA, ... )
}
