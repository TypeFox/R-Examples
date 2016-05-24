##' Show marker value sensitivity and specificity at mouse point
##'
##' A tag is set on the ROC curve at the mouse click and corresponding marker value, sensitivity and specificity shown below the click-point.
##' @title Click on ROC curve
##' @param object An object obtained with function \code{Roc}
##' @param pch the symbol of the tag
##' @param label If TRUE label the tag.
##' @param adj passed to \code{text} to adjust of the legend relative
##' to clickpoint.
##' @param col the color of the tag
##' @param cex the size of the tag
##' @param ... passed to \code{identify}
##' @return the values at the tag
##' @seealso identify Roc
##' @examples
##'
##'  \dontrun{
##' x <- abs(rnorm(20))
##' d <- data.frame(y=rbinom(1:20,1,p=x/max(x)))
##' r <- Roc(y~x,data=d)
##' plot(r)
##' click.Roc(r)
##'  }
##' @export 
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
click.Roc <- function(object, pch=19,label=TRUE, adj, col="orange", cex=3, ...) {
    x <- 1-object$Roc[[1]]$Specificity
    y <- object$Roc[[1]]$Sensitivity
    m <- object$breaks[[1]]
    n=length(x)
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, length(x)); res <- integer(0)
    while(sum(sel) < n) {
        ans <- identify(x[!sel], y[!sel], n=1, plot=FALSE, ...)
        if(!length(ans)) break
        ans <- which(!sel)[ans]
        spec <- 1-x[ans]
        sens <- y[ans]
        points(1-spec, sens, pch = pch, col=col, cex=cex)
        if (label) text(1-spec, sens, ans)
        if (missing(adj)){
            adj <- c(-.3,1)
        }
        text(1-spec,sens,adj=adj,paste(c("Marker=","Sensitivity=","Specificity="),c(round(m[ans],2),round(100*c(sens,spec),1)),collapse="\n"))
        sel[ans] <- TRUE
        res <- c(res, point=ans, marker=m[ans],Sensitivity=sens,Specificity=spec)
    }
    res
}
