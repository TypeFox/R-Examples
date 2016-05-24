#-------------------------------------------------
# This function has been taken from the package 
# "genalg" by E. Willighagen available on the CRAN
#-------------------------------------------------
plot.rbga <- function(x, type = "default", breaks = 10, ...) {
    rbga.object <- x
    if ((type == "trend") & (rbga.object$type == "floats chromosome")) {
        plot(rbga.object$best, type = "l", main = "", xlab = "Iteration (Generation)", 
            ylab = "Best (black lower line) and mean (red upper line) evaluation value", 
            ...)
        lines(rbga.object$mean, col = "red", ...)
    } else {
        if (type != "default") {
            warning(paste("Plot type", type, "not supported for a RBGA object of type", 
                rbga.object$type))
        }
        max <- max(rbga.object$best, rbga.object$mean)
        min <- min(rbga.object$best, rbga.object$mean)
        plot(rbga.object$best, type = "l", main = "", ylim = c(min, 
            max), xlab = "Iteration (Generation)", ylab = "Best (black lower line) and mean (red upper line) evaluation value", 
            ...)
        lines(rbga.object$mean, col = "red", ...)
    }
}
