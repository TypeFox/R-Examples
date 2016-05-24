library(adegraphics)
pdf("sarrow.pdf")

## ex1 : from tdr641
data(doubs, package = "ade4") 
dudi1 <- ade4::dudi.pca(doubs$env, scale = T, scan = F, nf = 3)
dudi2 <- ade4::dudi.pca(doubs$fish, scale = T, scan = F, nf = 2)
coin1 <- ade4::coinertia(dudi1, dudi2, scan = F, nf = 2)
g1 <- s.arrow(coin1$l1, plabels.cex = 0.87)
g2 <- s.arrow(coin1$c1, plabels.cex = 1)

## ex2 : from bs81
data(granulo, package = "ade4")
w <- data.frame(t(apply(granulo$tab, 1, function(x) x / sum(x))))
g3 <- s.arrow(ade4::dudi.pca(data.frame(w), scan = F, nf = 2)$co)

wtr <- data.frame(t(w))
wmoy <- data.frame(matrix(apply(wtr, 1, mean), 1))
dudi3 <- ade4::dudi.pca(w, scal = FALSE, scan = FALSE)
wmoy <- ade4::suprow(dudi3, wmoy)$lisup

g4 <- s.arrow(dudi3$c1, plabels.cex = 1.5)
g4 <- s.distri(dudi3$c1, wtr, starSize = 0.33, ellipseSize = 0, add = TRUE, plabels.cex = 1)
g4 <- s.label(wmoy, ppoint.cex = 5, plabels.cex = 0, add = TRUE)

## ex3
data(deug, package = "ade4")
pca1 <- ade4::dudi.pca(deug$tab, scal = FALSE, center = deug$cent, scan = FALSE)
g5 <- s.arrow(40 * pca1$c1)

## ex4
xy <- cbind(rnorm(50), rnorm(50))
g6 <- s.arrow(xy, plabels.cex = 0.9, parrows = list(angle = 20))
