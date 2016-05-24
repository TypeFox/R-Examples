barplot.noia.gpmap <-
function (height, GPcol = c("indianred", "palegreen", "royalblue"), 
    arrowscol = "purple", stderr = TRUE, main = NA, ylab = NA, 
    ...) 
{
    if (stderr) {
        fact <- 1
        err.name <- "stderr"
    }
    else {
        fact <- 1.96
        err.name <- "IC95%"
    }
    "barplotGPmap" <- function(matGval, err, ...) {
        lab1 <- c("aa", "aA", "AA")
        lab2 <- c("bb", "bB", "BB")
        mid <- barplot(matGval, col = GPcol, beside = TRUE, names.arg = lab1, 
            cex.main = 1, ylim = c(0, ceiling(max(matGval + err)) * 
                1.3), ...)
        if (is.matrix(matGval)) {
            legend("topleft", legend = lab2, pch = 15, col = GPcol, 
                bty = "n", cex = 0.8)
        }
        arrows(mid, matGval, mid, matGval + err, code = 2, length = 0.03, 
            angle = 90, col = arrowscol, cex.lab = 0.3)
    }
    col.lab <- c("cc", "cC", "CC", rep("", 6))
    y.lab <- c("dd", "", "", "dD", "", "", "DD", "", "")
    n <- nrow(height)/9
    if (n < 1) 
        n <- 1
    if (n > 3) 
        n <- 3
    par(mfrow = c(1, n), mar = c(2.2, 4, 1, 0.5))
    for (i in 1:n) {
        if (n == 3) {
            co <- col.lab[i]
            y <- "Genotypic values"
        }
        if (n > 3) {
            co <- col.lab[i]
            y <- y.lab[i]
        }
        else {
            co <- paste("GPmap", "(", err.name, ")")
            y <- "Genotypic values"
        }
        if (nrow(height) < 9) {
            matGval = height[, 1]
            err <- height[, 2] * fact
        }
        else {
            matGval <- matrix(height[((i * 9) - 8):(i * 9), 1], 
                3, 3)
            err <- (matrix(height[((i * 9) - 8):(i * 9), 2], 
                3, 3)) * fact
        }
        if (is.na(main)) {
            main <- co
        }
        if (is.na(ylab)) {
            ylab <- y
        }
        barplotGPmap(matGval, err, main = main, ylab = ylab, 
            ...)
    }
}
