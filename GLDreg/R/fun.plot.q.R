fun.plot.q <-
function (x, y, fit, quant.info, ...) 
{
    quant.info <-  quant.info[1:(nrow(quant.info)-2),,drop=F]
    ncol.quant.info <- ncol(quant.info)
    plot(x, y, ...)
    sapply(1:ncol.quant.info, function(i, quant.info, fit) {
        new.data <- cbind(x, data.matrix(fit$x) %*% quant.info[, 
            i])
        new.data <- new.data[order(new.data[, 1]), ]
        lines(new.data[, 1], new.data[, 2], col = "red")
    }, quant.info, fit)
}
