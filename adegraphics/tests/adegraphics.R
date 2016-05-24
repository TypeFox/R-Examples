library(adegraphics)
pdf("adegraphics.pdf")

xy <- cbind.data.frame(runif(7), runif(7))
g1 <- s.label(xy)

data(olympic, package = "ade4")
pca <- ade4::dudi.pca(olympic$tab, scan = FALSE)
g2 <- s.corcircle(pca$co, lab = names(olympic$tab))

g3 <- ADEgS(list(g1, g2), rbind(c(0, 0, 0.5, 1), c(0.5, 0, 1, 1)))
g4 <- ADEgS(list(g1, g2), layout = c(1, 2)) ## the same as g3
g4b <- ADEgS(list(g1, g2)) ## the same as g3
g5 <- s.label(xy, plabels.cex = 0, paxes.draw = TRUE, ppoints.col = "red")

g6 <- superpose(g1, g5, plot = TRUE)
g6b <- s.density(xy)
g7 <- superpose(s.density(xy), g5, plot = TRUE)
 
g8 <- superpose(s.label(xy, plabels.boxes.col = "orange", plot = FALSE), 
                s.label(xy, plabels.cex = 0, paxes.draw = TRUE, ppoints.col = "red", plot = FALSE), 
                plot = TRUE)

g9 <- g8[1, drop = TRUE]
class(g9)
g10 <- g8[1, drop = FALSE]
class(g10)

g11 <- ADEgS(list(g8, g3), positions = rbind(c(0, 0, 0.5, 1), c(0.5, 0, 1, 1)))

## cbindADEgS - rbindADEgS
g12 <- cbindADEg(g1, g2, plot = TRUE) ## the same as g3
g13 <- cbindADEg(g8, g3, plot = TRUE) ## the same as g11
g14 <- rbindADEg(g8, g3, plot = TRUE)

data(banque, package = "ade4")
banque.acm <- ade4::dudi.acm(banque, scann = FALSE, nf = 3)
g15 <- score(banque.acm, which = which(banque.acm$cr[, 1] > 0.2), plot = FALSE)
g15 <- g15[[1]]
cbindADEg(g15[[1]], g15[[2]], plot = TRUE)   ## work on trellis object

