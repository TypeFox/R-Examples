
library(gridGraphics)

require(grDevices)

matplot1 <- function() {
    matplot((-4:5)^2, main = "Quadratic") # almost identical to plot(*)
}

sines <- outer(1:20, 1:4, function(x, y) sin(x / 20 * pi * y))

matplot2 <- function() {
    matplot(sines, pch = 1:4, type = "o", col = rainbow(ncol(sines)))
}

matplot3 <- function() {
    matplot(sines, type = "b", pch = 21:23, col = 2:5, bg = 2:5,
            main = "matplot(...., pch = 21:23, bg = 2:5)")
}

x <- 0:50/50

matplot4 <- function() {
    matplot(x, outer(x, 1:8, function(x, k) sin(k*pi * x)),
            ylim = c(-2,2), type = "plobcsSh",
            main= "matplot(,type = \"plobcsSh\" )")
}

matplot5 <- function() {
    ## pch & type =  vector of 1-chars :
    matplot(x, outer(x, 1:4, function(x, k) sin(k*pi * x)),
            pch = letters[1:4], type = c("b","p","o"))
}

matplot6 <- function() {
    lends <- c("round","butt","square")
    matplot(matrix(1:12, 4), type="c", lty=1, lwd=10, lend=lends)
    text(cbind(2.5, 2*c(1,3,5)-.4), lends, col= 1:3, cex = 1.5)
}

table(iris$Species) # is data.frame with 'Species' factor
iS <- iris$Species == "setosa"
iV <- iris$Species == "versicolor"

matplot7 <- function() {
    par(bg = "bisque")
    matplot(c(1, 8), c(0, 4.5), type =  "n", xlab = "Length", ylab = "Width",
            main = "Petal and Sepal Dimensions in Iris Blossoms")
    matpoints(iris[iS,c(1,3)], iris[iS,c(2,4)], pch = "sS", col = c(2,4))
    matpoints(iris[iV,c(1,3)], iris[iV,c(2,4)], pch = "vV", col = c(2,4))
    legend(1, 4, c("    Setosa Petals", "    Setosa Sepals",
                   "Versicolor Petals", "Versicolor Sepals"),
           pch = "sSvV", col = rep(c(2,4), 2))
}

nam.var <- colnames(iris)[-5]
nam.spec <- as.character(iris[1+50*0:2, "Species"])
iris.S <- array(NA, dim = c(50,4,3),
                dimnames = list(NULL, nam.var, nam.spec))
for(i in 1:3) iris.S[,,i] <- data.matrix(iris[1:50+50*(i-1), -5])

matplot8 <- function() {
    matplot(iris.S[, "Petal.Length",], iris.S[, "Petal.Width",], pch = "SCV",
            col = rainbow(3, start = 0.8, end = 0.1),
            sub = paste(c("S", "C", "V"), dimnames(iris.S)[[3]],
                sep = "=", collapse= ",  "),
            main = "Fisher's Iris Data")
}

plotdiff(expression(matplot1()), "matplot-1")
plotdiff(expression(matplot2()), "matplot-2")
plotdiff(expression(matplot3()), "matplot-3")
plotdiff(expression(matplot4()), "matplot-4")
plotdiff(expression(matplot5()), "matplot-5")
plotdiff(expression(matplot6()), "matplot-6")
plotdiff(expression(matplot7()), "matplot-7")
plotdiff(expression(matplot8()), "matplot-8")

plotdiffResult()
