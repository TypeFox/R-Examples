"plot.drift" <-
function(x, main = NULL, xlab = NULL, ylab = NULL, ...){
    if (!((inherits(x, "bounds"))|(inherits(x, "drift"))))
     stop("'x' must inherit from class \"bounds\" or \"drift\"")
    if (is.null(main))
       main <- "Sequential boundaries using the Lan-DeMets method"
    if (is.null(xlab))
       xlab <- "Time"
    z <- x$time
    r <- array(0, dim=length(z))
    if (inherits(x, "bounds"))
       if (x$bounds.type==1){
          u <- x$upper.bounds
          z <- x$time
          ans <- xyplot(u+r~z, main = main, xlab = xlab, ylab = ylab, type=c("p","l"), col = "black", ...)
       }
       else{
          u <- x$upper.bounds
          l <- x$lower.bounds
          z <- x$time
          ans <- xyplot(u+l+r~z, main = main, xlab = xlab, ylab = ylab, type=c("p","l"), col = "black", ...)
       }
    if  (inherits(x, "drift")){
       u <- x$upper.bounds
       l <- x$lower.bounds
       ans <- xyplot(u+l+r~z, main = main, xlab = xlab, ylab = ylab, type=c("p","l"), col = "black", ...)
    }
    ans
}

