#' Plot Jn values
#'
#' Plot the \eqn{J_n(x)} values output from \code{\link{snowball}}, with significant genes highlighted. See references for more details.
#' @param x an output from \code{\link{snowball}}
#' @param fs the corresponding output from \code{\link{select.features}}
#' @param pch.nonsig \code{pch} of the symbols for non-significant genes. See \code{\link{par}} for more details
#' @param pch.sig \code{pch} of the symbols for significant genes.
#' @param below.median a \code{logical} value, set to \code{TRUE} if the genes blow the median are to be highlighted 
#' @param col.above set the highlight color for genes above the median
#' @param col.below set the highlight color for genes below the median
#' @export
plotJn <- function(x,fs, pch.nonsig=21,pch.sig=19,
		   below.median=T,
		   col.above="red",col.below="red")
{
    plot(seq(along=x$weights),
         x$weights,
         xlab="Gene Index",
         ylab="Jn score",
         main="Snowball analysis",
         type="n")
    sigs <- row.names(x)%in%row.names(fs$selectedList)
    positives <-
      row.names(x) %in% row.names(subset(fs$selectedList,subset=fs$selectedList$positive))
    points(seq(along=x$weights)[!sigs],x$weights[!sigs],pch=pch.nonsig)
    points(seq(along=x$weights)[positives],x$weights[positives],pch=pch.sig,col=col.above)
    negatives <-	
	row.names(x) %in% row.names(subset(fs$selectedList,subset=!fs$selectedList$positive))

    if(below.median) {
	points(seq(along=x$weights)[negatives],x$weights[negatives],pch=pch.sig,col=col.below)
    } else {
	points(seq(along=x$weights)[negatives],x$weights[negatives],pch=pch.nonsig)
    }
}
