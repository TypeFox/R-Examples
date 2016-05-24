#' Area vs. x-axis
#'
#' Function charts density plot of a single spot following x-axis.
#' @param object S3 object of working TLC
#' @param spot Number of the spot (counted left to right).
#' @param plot Boolean, TRUE default and displays densitometric distribution.
#' @param returndf Boolean, TRUE by default, returns \code{data.frame} with \code{x} and \code{Area} values.
#' @param ptype Point type for the plot. Default "o" (Uses same values as \code{type} variable from \code{plot} function)
#' @param ... Additional parameters (for \code{plot} type function).
#'
#' @return Returns \code{data.frame} with \code{x} and \code{Area} values.
#'
#' @author Ivan D. Pavicevic, \email{ivanp84@@gmail.com}
#' 
#' @examples
#' 
#' # A test example
#' fname01 <- system.file("extdata", "test025to100sp.tiff", package="qtlc")
#' testTLC <- createTLC(fname01, RGB=FALSE)
#' print(testTLC)
#'
#' # now we'll imitate interactive spot2D function,
#' # and create spots coordinates automatically,
#' # for interactive version run testTLC <- spot2D(testTLC)
#' testTLC$spots$x <- c(40.93354, 83.18687, 121.59899, 160.01111, 203.54485,
#'                      239.39616, 280.36909, 320.06161, 362.31494, 399.44666,
#'                      439.13919, 480.11211, 518.52423, 559.49716, 599.18969)
#' testTLC$spots$y <- c(198.3160, 198.3160, 199.2833, 198.3160, 198.3160,
#'                      198.3160, 198.3160, 198.3160, 197.3487, 198.3160,
#'                      199.2833, 198.3160, 199.2833, 199.2833, 199.2833)
#'
#' testTLC <- select2D(testTLC, 30, 30)
#' testTLC <- matrices2D(testTLC)
#' testTLC <- summat2D(testTLC)
#'
#' # and now test the areadens2D for each spot
#' par(mfrow=c(3,3))
#' for(i in 1:15) {
#' areadens2D(testTLC, spot=i, ptype="l")
#' }
#'
#' @export
#' @importFrom graphics plot
#' 
areadens2D <- function(object, spot=NULL, plot=TRUE, returndf = TRUE, ptype="o", ...) {
    if (is.null(spot)) {
        cat("\nSpecify the spot for density plot.\n");
        return(NULL);
    }

    mat <- object$spot_matrices[,,spot];
    P <- NULL;
    
    for(i in 1:nrow(mat)) {
        P <- c(P, sum(mat[i,]));
    }

    Pdf <- data.frame(1:nrow(mat), P);
    names(Pdf) <- c("x", "Area");
    if (plot == TRUE) {graphics::plot(Pdf, type=ptype, main=spot, ...);}
    if (returndf == TRUE) {return(Pdf);}
}

