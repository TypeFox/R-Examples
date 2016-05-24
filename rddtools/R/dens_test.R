#' Run the McCracy test for manipulation of the forcing variable
#' 
#' Calls the \code{\link[rdd]{DCdensity}} test from package \code{rdd} on a \code{rdd_object}.
#' 
#' @param rdd_object object of class rdd_data
#' @param bin Argument of the \code{\link{DCdensity}} function, the binwidth
#' @param bw Argument of the \code{\link{DCdensity}} function, the bandwidth
#' @param plot Whether to return a plot. Logical, default ot TRUE. 
#' @param \ldots Further arguments passed to \code{\link[rdd]{DCdensity}}. 
#' @export
#' @import rdd
#' @examples
#' data(house)
#' house_rdd <- rdd_data(y=house$y, x=house$x, cutpoint=0)
#' dens_test(house_rdd)



dens_test <- function(rdd_object, bin = NULL, bw = NULL, plot = TRUE, ...) {
    checkIsRDD(rdd_object)
    cutpoint <- getCutpoint(rdd_object)
    x <- getOriginalX(rdd_object)
    test <- try(DCdensity(runvar = x, cutpoint = cutpoint, bin = bin, bw = bw, plot = plot, ext.out = TRUE, ...), silent = TRUE)
    if (inherits(test, "try-error")) {
        warning("Error in computing the density, returning a simple histogram", if (is.null(bin)) 
            " with arbitrary bin" else NULL)
        if (is.null(bin)) {
            test <- try(DCdensity(rdd_object$x, cutpoint, bin = bin, bw = 0.2, ext.out = TRUE, plot = FALSE), silent = TRUE)
            bin <- test$binsize
        }
        max_x <- max(rdd_object$x, na.rm = TRUE)
        seq_breaks <- seq(from = min(rdd_object$x, na.rm = TRUE), to = max_x, by = bin)
        if (max_x > max(seq_breaks)) 
            seq_breaks <- c(seq_breaks, max_x + 0.001)
        hist(rdd_object$x, breaks = seq_breaks)
        abline(v = cutpoint, col = 2, lty = 2)
    }
    
    test.htest <- list()
    test.htest$statistic <- c(`z-val` = test$z)
    test.htest$p.value <- test$p
    test.htest$data.name <- deparse(substitute(rdd_object))
    test.htest$method <- "McCrary Test for no discontinuity of density around cutpoint"
    test.htest$alternative <- "Density is discontinuous around cutpoint"
    test.htest$estimate <- c(Discontinuity = test$theta)
    test.htest$test.output <- test
    class(test.htest) <- "htest"
    return(test.htest)
} 
