"summary.bounds" <-
function (object, ...) 
{
    z <- object
    if (!inherits(z, "bounds")) 
     stop("'object' must inherit from class \"bounds\"")
    p <- length(z$time)
    if (identical(z$time,z$time2)){
       b <- matrix(NA, p, 5)
       b[,1:5] <- c(z$time, z$lower.bounds, z$upper.bounds, z$exit.pr, z$diff.pr) 
       colnames(b) <- c("Time", "Lower", "Upper", "Exit pr.", "Diff. pr.")
    }
    else{
       b <- matrix(NA, p, 6)
       b[,1:6] <- c(z$time, z$time2, z$lower.bounds, z$upper.bounds, z$exit.pr, z$diff.pr) 
       colnames(b) <- c("Time", "Time 2", "Lower", "Upper", "Exit pr.", "Diff. pr.")
    }
    ans <- list()
    ans$type <- z$bounds.type
    ans$spending <- z$spending.type
    ans$n <- p
    ans$alpha <- z$alpha    
    ans$oalpha <- z$overall.alpha
    ans$bounds <- b
    class(ans) <- "summary.bounds"
    return(ans)
}

