`performance.wa` <- function(object, ...) {
    retval <- with(object, c(rmse, r.squared, avg.bias,
                             max.bias))
    names(retval) <- c("RMSE","R2","Avg.Bias","Max.Bias")
    class(retval) <- "performance"
    retval
}
