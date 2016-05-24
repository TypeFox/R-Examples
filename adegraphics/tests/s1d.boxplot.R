library(adegraphics)
pdf("s1d.boxplot.pdf")

## ex1
x <- c(rnorm(10), rnorm(10))
fact <- factor(rep(c("A", "B"), 10))
g1 <- s1d.boxplot(x, fact)
g2 <- s1d.boxplot(x, fact, ppolygon.border = c("red", "blue"), box.rectangle = list(alpha = 1, fill = "green"))

## ex2
w1 <- rnorm(100, -1)
w2 <- rnorm(100)
w3 <- rnorm(100, 1)
f1 <- gl(3, 100)
f2 <- gl(30, 10)
g3 <- s1d.boxplot(c(w1, w2, w3), f1)
g4 <- s1d.boxplot(c(w1, w2, w3), f2)
g5 <- s1d.boxplot(c(w1, w2, w3), f2, p1d.rug.draw = FALSE)

mat <- matrix(0, ncol = 1, nrow = 8)
mat[c(2), ] <- 1
mat[c(3:8), ] <- 2
mat[1, ] <- 3
g6 <- ADEgS(c(g3, g4, s1d.label(c(w1, w2, w3), p1d = list(rug = list(tck = 0.8), rev = TRUE), ppoints.cex = 0, plabels.cex = 0, plot = F, pgrid.draw = F)), 
  layout = matrix((rev(mat)), ncol = 1))
g7 <- s1d.boxplot(c(w1, w2, w3), data.frame(f1, f2))

## ex3
data(banque, package = "ade4")
banque.acm <- ade4::dudi.acm(banque, scan = FALSE, nf = 4)
s1d.boxplot(banque.acm$l1[, 1], banque[, 1:7], plabels.cex = 1.8)
