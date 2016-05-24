plot.demosimhybrid <-
structure(function (x, col = c(2, 3, 4, "orange", "orchid", 7),...) 
{
    l <- dim(x)[1]
    g <- rep("G", l)
    num <- (0:(l - 1))
    gen <- paste(g, num)
    layout(matrix(c(1, 2), byrow = F), heights = c(0.1, 0.5))
    par(mar = c(1, 1, 1, 1))
    plot(1, type = "n", axes = FALSE, ann = FALSE, col = c(2, 
        3, 4, 0, 7, 4))
    legend(x = "center", legend = c("PA", "PB", expression("F"[1]), "BxA", 
        "BxB", "Fx"), fill = col, cex = 1.5, bty = "n", xjust = 0, 
        yjust = 0, ncol = 5)
    par(mar = c(4, 6, 1, 2))
    barplot(t(x), names.arg = gen, col = col, ylab = "Proportion of individuals", 
        cex.lab = 1.4)
}, export = FALSE, S3class = "demosimhybrid", modifiers = "public")
