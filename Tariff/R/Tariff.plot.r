#' Plot CSMF of the results obtained from Tariff algorithm
#'
#' This function plots the CSMF of the fitted results.
#'
#' @param x fitted object from \code{\link{tariff}}
#' @param top maximum causes to plot
#' @param min.prob minimum fraction for the causes plotted
#' @param ... Arguments to be passed to/from graphic function
#'
#' @examples
#'
#'\donttest{
#' data("RandomVA3")
#' test <- RandomVA3[1:200, ]
#' train <- RandomVA3[201:400, ]
#' allcauses <- unique(train$cause)
#' fit <- tariff(causes.train = "cause", symps.train = train, 
#' 				symps.test = test, causes.table = allcauses)
#' plot(fit, top = 10, main = "Top 5 population COD distribution")
#' plot(fit, min.prob = 0.05, main = "Ppulation COD distribution (at least 5%)")
#' }

plot.tariff<- function(x, top = NULL, min.prob = 0, ...){
	
	sx <- summary(x)
	dist.cod <- sx$csmf 

    if(!is.null(top)){
        if(top < length(dist.cod)){
            thre <- sort(dist.cod, decreasing=TRUE)[top]
            min.prob <- max(min.prob, thre)
        }
    }
    dist.cod.min <- dist.cod[dist.cod >= min.prob ]
    dist.cod.min <- sort(dist.cod.min, decreasing = FALSE)
    par(las = 2)
    par(mar = c(5,15,4,2))
    bar.color <- grey.colors(length(dist.cod.min))
    bar.color <- rev(bar.color)
    barplot(dist.cod.min , horiz = TRUE,names.arg = names(dist.cod.min), col = bar.color, cex.names=0.8, xlab = "Probability", ...)

}