"fun.plot.fit.bm" <-
function (fit.obj, data, nclass = 50, xlab = "", name = "", main = "", 
    param.vec, ylab = "Density") 
{
    fit.obj <- fit.obj$par
    prop <- fit.obj[9]
    hist.lim <- histsu(data, probability = TRUE, plot = FALSE, nclass = nclass)$density
    first.part <- prop * dgl(seq(min(data), max(data), length = 1000), 
        fit.obj[1], fit.obj[2], fit.obj[3], fit.obj[4], param = param.vec[1])
    second.part <- (1 - prop) * dgl(seq(min(data), max(data), 
        length = 1000), fit.obj[5], fit.obj[6], fit.obj[7], fit.obj[8], 
        param = param.vec[2])
    jmax <- max(hist.lim, first.part + second.part)
    if (main != "") {
        histsu(data, probability = TRUE, xlab = xlab, ylab = ylab, 
            nclass = nclass, ylim = c(0, jmax), main = main)
        lines(seq(min(data), max(data), length = 1000), first.part, 
            col = 5, lwd = 4)
        lines(seq(min(data), max(data), length = 1000), second.part, 
            col = 5, lwd = 4)
        lines(seq(min(data), max(data), length = 1000), first.part + 
            second.part, col = 4, lwd = 4)
        if (name != "") {
            legend("topright", name, lwd = 4, col = 5)
        }
    }
    if (main == "") {
        histsu(data, probability = TRUE, xlab = xlab, ylab = ylab, 
            nclass = nclass, ylim = c(0, jmax), main = paste(name, 
                param.vec[1], param.vec[2]))
        lines(seq(min(data), max(data), length = 1000), first.part, 
            col = 5, lwd = 4)
        lines(seq(min(data), max(data), length = 1000), second.part, 
            col = 5, lwd = 4)
        lines(seq(min(data), max(data), length = 1000), first.part + 
            second.part, col = 4, lwd = 4)
    }
}

