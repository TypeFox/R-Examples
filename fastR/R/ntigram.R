#' @export
ntigram <-
function (x, bins = 12, name = deparse(substitute(x)), ...) 
{
    if (bins >= length(x)/3) {
        bins <- round(length(x)/3)
    }
    sdata <- sort(x)
    width <- round(length(sdata)/bins)
    l <- length(sdata)
    cdata <- sdata[seq(round(l/(2 * bins)), l, length = bins)]
    bks <- cdata[0]
    for (j in 2:length(cdata)) {
        if (cdata[[j]] != cdata[[j - 1]]) {
            bks <- c(bks, cdata[j])
        }
    }
    top <- 2 * max(x) - max(x[x < max(x)])
    bot <- 2 * min(x) - min(x[x > min(x)])
    bks <- c(bot, bks, top)
    cat("break pts for ntigram = ", bks, "\n")
    histogram(x, breaks = bks, main = paste("Histogram of", name), 
        xlab = name, type = "density", ...)
}
