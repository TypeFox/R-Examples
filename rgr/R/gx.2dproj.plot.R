gx.2dproj.plot <-
function (save, rowids = NULL, main = "", pch = 3, cex = 0.8, col = 1, ...) 
{
    frame()
    if (main == "") 
        banner <- save$main
    else banner <- main
    if(is.null(rowids)) {
        plot(save$x, save$y, xlab = save$xlab, ylab = save$ylab, 
            main = banner, pch = pch, cex = cex, col = col, ...)
    }
    else {
        plot(save$x, save$y, xlab = save$xlab, ylab = save$ylab, 
            type = "n", main = banner, ...)
        if(rowids) text(save$x, save$y, save$row.numbers, cex = cex, col = col, ...)    
        else text(save$x, save$y, save$matnames[[1]], cex = cex, col = col, ...)
    }
    abline(v = 0, lty = 2)
    abline(h = 0, lty = 2)
    invisible()
}
