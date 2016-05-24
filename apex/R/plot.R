
######################
####  PLOT METHOD ####
######################
#' Display multidna objects
#'
#' Default printing for multidna objects
#'
#' @export
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @param x a multidna object
#' @param y an integer vector indicating the genes to plot
#' @param rows a logical indicating if different genes should be displayed in separate rows
#' @param ask a logical indicating if the user should be prompted between graphs
#' @param ... arguments passed to \code{\link{image.DNAbin}}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @import ape
#'
#' @aliases plot,multidna-method
#' @aliases plot.multidna
#'
#' @docType methods
#'
#' @examples
#' ## simple conversion with nicely ordered output
#' data(woodmouse)
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
#' x <- new("multidna", genes)
#' x
#' plot(x)
#'
setMethod ("plot", "multidna", function(x, y, rows=TRUE, ask=FALSE, ...){
    ## HANDLE ARGUMENTS ##
    n.genes <- length(x@dna)
    if(missing(y)) y <- 1:n.genes
    y <- as.integer(y)
    y <- y[y>0 | y<=n.genes]

    ## MAKE PLOT ##
    opar <- par(no.readonly=TRUE)
    on.exit(par(opar))
    par(ask=ask)
    if(rows) par(mfrow=c(n.genes,1))
    for(i in y){
        image(x@dna[[i]], ...)
        mtext(side=1, text=names(x@dna)[i], line=3, cex=2)
    }
})





