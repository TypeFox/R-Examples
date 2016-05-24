library(adegraphics)
pdf("s1d.gauss.pdf")

data(meau, package= "ade4")
envpca <- ade4::dudi.pca(meau$env, scannf = FALSE)
dffac <- cbind.data.frame(meau$design$season, meau$design$site)
g1 <- s1d.gauss(envpca$li[, 1], dffac[, 1])
g2 <- s1d.gauss(envpca$li[, 1], dffac[, 1], ppoly.col = 1:4, fill = TRUE, plines.col = 1:4, col = FALSE)

g3 <- s1d.gauss(envpca$li[, 1], dffac[, 2], ppoly.col = 1:4, paxes.draw = TRUE, ylim = c(0, 2), fill = T, p1d.hori = F)
update(g3, p1d.rev = T)

g4 <- s1d.gauss(envpca$li[, 1], fac = dffac, fill = TRUE, col = 1:5)
g5 <- s1d.gauss(envpca$li[, 1], fac = dffac, fill = T, col = F, ppoly.col = 1:6)
g6 <- s1d.gauss(envpca$li[, 1], fac = dffac[, 1], fill = T, col = 1:6, ppoly.col = 1:6)

g7 <- s1d.gauss(envpca$li[, 1], fac = dffac, fill = T, col = 1:6, ppoly.col = 1:6, steps = 10)
g8 <- s1d.gauss(envpca$li[, 1], dffac[, 2])
