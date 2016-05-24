library(adegraphics)
pdf("smatch.pdf")

X <- data.frame(x = runif(50, -1, 2), y = runif(50, -1, 2))
Y <- X + rnorm(100, sd = 0.3)

g1 <- s.match(X, Y, ppoints.cex = 0, col = c("blue", "red"))
g2 <- s.match(X, Y, arr = FALSE, ppoints.cex = 2, ppoints.col = c("blue", "green"))
g3 <- s.match(X, Y, arr = FALSE)
g4 <- s.match(X, Y, arrows = TRUE, plabels = list(alpha = 1, col = "black", cex = 1), plines = list(col = "red"), panel.background = list(col = "antiquewhite"))
