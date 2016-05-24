library(adegraphics)
pdf("s1d.distri.pdf")

## ex1
score <- rnorm(10)
df <- data.frame(matrix(rep(c(1, 2), 10), ncol = 2))
g1 <- s1d.distri(score, df)

## ex2
set.seed(1)
w <- seq(-1, 1, le = 200)
distri <- data.frame(lapply(1:50, function(x) sample(200:1) * ((w >= (- x / 50)) & (w <= x / 50))))
names(distri) <- paste("w", 1:50, sep = "")
g2 <- s1d.distri(w, distri)
g3 <- s1d.distri(w, distri, yrank = TRUE, sdS = 1.5)
