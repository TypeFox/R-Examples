#' @name autoplotTableAndPlot
#' @title Arrange a survival plot with corresponding table and legend.
#' 
#' @include autoplotTen.R
#' 
#' @param object An object of class \code{"tableAndPlot"}, as returned by
#' \code{ggplot.Ten}.
#' @param ... Additional arguments (not implemented).
#' @param hideTabLeg Hide table legend.
#' \cr
#' If \code{hideTabLeg = TRUE} (the default), the table legend will not appear.
#' @param tabHeight Table height, as a fraction/ proportion of the whole.
#'  \cr
#' \code{tabHeight=0.25} (the default) makes the table 
#' \eqn{0.25 = 25\%} of the whole plot height.
#' 
#' @return A graph, plotted with \code{gridExtra::grid.arrange}.
#' 
#' @details Arguments to \code{plotHeigth} and \code{tabHeight} are
#' best specified as fractions adding to \eqn{1},
#' \cr
#' 
#' @note This method is called by \code{\link{print.tableAndPlot}}
#' and by \code{print.stratTableAndPlot}.
#' 
#' @author Chris Dardis. Based on existing work by
#' R. Saccilotto, Abhijit Dasgupta, Gil Tomas and Mark Cowley.
#'
#' @keywords graphics
#'
#' @rdname autoplotTAP
#' @method autoplot tableAndPlot
#' @aliases autoplot.tableAndPlot
#' @export
#' @examples
#' data("kidney", package="KMsurv")
#' autoplot(survfit(Surv(time, delta) ~ type, data=kidney), type="fill")
#' autoplot(ten(survfit(Surv(time, delta) ~ type, data=kidney)), type="fill")
#' data("bmt", package="KMsurv")
#' s2 <- survfit(Surv(time=t2, event=d3) ~ group, data=bmt)
#' autoplot(s2)
#' 
autoplot.tableAndPlot <- function(object, ...,
                                  hideTabLeg=TRUE,
                                  tabHeight=0.25){
    stopifnot(inherits(object, "tableAndPlot"))
    stopifnot(0 < tabHeight & tabHeight < 1)
    if (hideTabLeg) {
        object$table <- object$table +
            theme(legend.key.height=NULL,
                  legend.key.width=NULL,
                  legend.key=element_rect(colour=NA, fill=NA),
                  legend.text=element_text(colour=NA),
                  legend.title=element_text(colour=NA))
    }
    ## change to graphical objects
    grobs1 <- lapply(rev(seq.int(object)),
                     function(i) ggplotGrob(object[[i]]))
    ## collect the widths for each grob of each plot
    w1 <- lapply(seq.int(grobs1),
                 function(i) grobs1[[i]]$widths[2:5])
    ## use do.call to get the max width
    maxWidth1 <- do.call(grid::unit.pmax, w1)
    ## asign the max width to each grob
    for (i in seq.int(grobs1)) {
        grobs1[[i]]$widths[2:5] <- as.list(maxWidth1)
    }
    ## plot
    do.call(gridExtra::grid.arrange, c(grobs1,
                                       nrow=2,
                                       heights=list(c(1 - tabHeight, tabHeight))))
}
