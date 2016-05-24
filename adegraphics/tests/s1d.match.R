library(adegraphics)
pdf("s1d.match.pdf")

g1 <- s1d.match(-5:5, 2 * (-5:5))
g2 <- s1d.match(rnorm(10), runif(10), p1d.hor = FALSE)
g3 <- s1d.match(1:5, 7:11, p1d.hor = F, p1d.rev = T)
