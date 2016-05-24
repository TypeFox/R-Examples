################################################################################
# Random KNN Plot Functions                                                    #
# File:   plot.R                                                               #
# Author: Shengqiao Li                                                         #
# Date:   June 24, 2008 (initial)                                              #
################################################################################
plot.rknnBeg <-function (x, col = "springgreen4", xlab = "no. of features", ylab = "mean accuracy", ...)
{
    if (!inherits(x, "rknnBeg"))
         stop(deparse(substitute(x)), " is not a rknnBeg object")
         
    plot(x$p, x$mean_accuracy, xlim = rev(range(x$p)), type = "s", 
        col = col, log = "x", xlab =xlab, ylab = ylab, ...);
}
plot.rknnBel <- function (x, col = "springgreen4", xlab = "no. of features", ylab = "mean accuracy", ...)
{
    if (!inherits(x, "rknnBel")) stop(deparse(substitute(x)), " is not a rknnBel object")
    plot(x$p, x$mean_accuracy, xlim = rev(range(x$p)), type = "s", col = col, xlab = xlab, ylab = ylab, ...)
}
plot.rknnSupport <- function (x, n.var = min(30, length(x$support)), 
                  main = deparse(substitute(x)), bg="gold", lcolor="blue",  ...)
{
    if (!inherits(x, "rknnSupport")) 
        stop(deparse(substitute(x)), " is not a rknnSupport object")
        
    dotchart(rev(x$support[1:n.var]), xlim = rev(range(x$support[1:n.var])), 
        xlab = "Support", bg = bg, lcolor = lcolor, main = main, ...);
}
################################################################################
