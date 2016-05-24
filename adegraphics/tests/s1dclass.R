library(adegraphics)
pdf("s1dclass.pdf")
data(meau, package = "ade4")

## ex1
envpca <- ade4::dudi.pca(meau$env, scannf = FALSE)
g1 <- s1d.class(envpca$li[, 1], poslab = "value", meau$design$season, col = 1:6)
update(g1, p1d.horizontal = FALSE)
update(g1, p1d.reverse = TRUE)
g2 <- s1d.class(envpca$li[, 1], meau$design$season, col = 1:6, p1d.reverse = TRUE)
g3 <- s1d.class(envpca$li[, 1], meau$design$season, col = 1:6, p1d.hori = F)

## ex2
set.seed(0)
score1 <- c(rnorm(3, mean = 0, sd = 0.5), rnorm(3, mean = 1, sd = 0.5), rnorm(5, mean = 2, sd = 0.5))
factor1 <- factor(rep(LETTERS[1:3], times = c(3, 3, 5)))
g4 <- s1d.class(score1, factor1, col = 1:3)

## ex3
score2 <- c(rnorm(10, mean = 0, sd = 0.5), rnorm(15, mean = -1, sd = 0.2), rnorm(10, mean = 2, sd = 0.5))
factor2 <- factor(rep(c(1, 3, 2), times = c(10, 15, 10)))
levels(factor2) <- c("mean0", "mean2", "mean-1")
g5 <- s1d.class(score2, factor2, col = 1:3)
update(g5, posla = "value")

indx <- rank(rnorm(35))
factor2 <- factor2[rank(indx)]
s1d.class(score2[indx], factor2[indx], col = 1:3, posla = "regular")
