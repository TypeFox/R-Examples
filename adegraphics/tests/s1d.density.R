library(adegraphics)
pdf("s1d.density.pdf")

x1 <- rnorm(1000)
g1 <- s1d.density(x1)
g2 <- s1d.density(x1, fill = T, ppoly.col = "blue", p1d.rev = T)

x2 <- c(rnorm(1000, mean = -0.5, sd = 0.5), rnorm(1000, mean = 1))
fact <- rep(c("A", "B"), each = 1000)
g3 <- s1d.density(x2, fact, col = c("red", "blue"))
g4 <- s1d.density(x2, fact, col = FALSE, ppoly.col = 2:3, fill = T, p1d.rev = F)
g5 <- s1d.density(x2, fact, col = FALSE, ppoly.col = 2:3, p1d.horizontal = F, p1d.rev = F, fill = TRUE)
