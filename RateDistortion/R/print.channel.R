print.channel <-
function(x, ...) {
    ## S3 method for displaying a channel

    channel <- x
    cat("\nCapacity limited information channel.\n")

    nj <- nrow(channel$x)
    dim.x <- ncol(channel$x)

    nk <- nrow(channel$y)
    dim.y <- ncol(channel$y)
    
    cat("  Dimension = [", nj, " x ", dim.x,
        "] --> [", nk, " x ", dim.y, "]\n", sep = "")

    if(abs(channel$R - 1) < 1e-10) {
        plural <- ""
    } else {
        plural <- "s"
    }
    
    cat("  Rate = ", channel$R, " bit",
        plural, "\n  Distortion = ", channel$D, "\n", sep = "")
    if(channel$termination == 0) {
        cat("\nBlahut algorithm converged to eps =",
            channel$eps, "in", channel$iters, "iterations.\n\n")
    } else if(channel$termination == 1) {
        cat("Algorithm exceeded maximum number of iterations.\n\n")
    }
}
