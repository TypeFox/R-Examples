`plot.snp` <-
function (x, type = barplot, label, ...)
{
    if (!inherits(x, "snp")) 
        stop("snp must be an object of class 'WGassociation'")
    if (missing(label)) 
        label <- deparse(substitute(x))
    old.mar <- par("mar")
    old.mfrow <- par("mfrow")
    on.exit(par(mar = old.mar, mfrow = old.mfrow))
    m <- m <- matrix(c(1, 2), nrow = 2, ncol = 1, byrow = TRUE)
    layout(m, heights = c(1, 5.5))
    par(mar = c(0, 0, 0, 0))
    xx <- summary(x)
    plot(c(1:5), rep(1, 5), ylim=c(0.1,1.6), type = "n", axes = FALSE, xlab = "", 
        ylab = "")
    text(1, 1.5, label, font = 2, adj = 0)
    crea.lab(xx$allele.freq, 1.6, 0.8, 0.25)
    crea.lab(xx$genotype.freq, 2.8, 0.8, 0.25)
    text(4.5, 1, paste("HWE (pvalue):", round(xx$HWE, 6)), cex = 0.8)
    par(mar = old.mar)
    type(xx$genotype.freq[, 1], ...)
}

