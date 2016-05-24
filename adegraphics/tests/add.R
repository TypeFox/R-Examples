library(adegraphics)
pdf("add.pdf")

data(granulo, package = "ade4")
df <- data.frame(t(apply(granulo$tab, 1, function(x) x / sum(x))))
pca <- ade4::dudi.pca(df, scal = FALSE, scan = FALSE)

g1 <- s.arrow(ade4::dudi.pca(data.frame(df), scan = F, nf = 2)$co)
g2 <- s.label(pca$li, plabels.cex = 0.5, plabels.col = "blue", plot = F)
g3 <- add.ADEg(g2)
g4 <- s.label(pca$c1, plabels.col = "red", add = T)
g5 <- s.arrow(pca$c1, plabels.cex = 1.5)
g6 <- ADEgS(list(g1 = g1, g5 = g5), layout = c(1, 2))
update(g6, pback.col = "lightblue", g1.plabels.cex = 2, g5.plabels.col = "red")


## addtext
# on a ADEg
addtext(g1, -1, 1, "Data Granulo", plabels.cex = 1.5, plabels.col = "red")
addtext(g1, -1, 1, c("Data", "Granulo"), plabels.cex = 1.5, plabels.col = c("red", "blue"))
addtext(g1, -1, 1, c("Data", "Granulo"), plabels.cex = 1.5, plabels.col = "red")

addtext(g1, c(-1, -0.5), 1, "Data Granulo", plabels.cex = 1.5, plabels.col = c("red", "blue"))
addtext(g1, c(-1, -0.5), 1, c("Data", "Granulo"), plabels.cex = 1.5, plabels.col = c("red", "blue"))

addtext(g1, -1, c(1, 0.9), c("Data", "Granulo"), plabels.cex = 1.5, plabels.col = "red")
addtext(g1, -1, c(1, 0.9), c("Data", "Granulo"), plabels.cex = c(1.5, 2), plabels.col = "red")

data(dunedata, package = "ade4")
afc1 <- ade4::dudi.coa(dunedata$veg, scannf = FALSE)
g7 <- table.value(dunedata$veg, symbol = "circle", ppoints.cex = 0.5, plot = FALSE)
addtext(g7, 1, 20, "A", plabels.srt = 45)

# on a trellis
library(lattice)
g8 <- xyplot(1 ~ 1)
addtext(g8, c(1.1, 0.9), 1.1, "trellis test", plabels.col = c("orange", "grey"))

# on a ADEgS
addtext(g6, -0.9, -0.1, "Data Granulo", plabels.cex = 1.2, which = 1)
addtext(g6, -0.5, -0.1, "Data Granulo", plabels.cex = 1.2, which = 2)
addtext(g6, -0.9, -0.1, "Data Granulo", plabels.cex = 1.2, which = 1:2)
addtext(g6, c(-0.9, -0.5), c(-0.2, -0.1), "Data Granulo", plabels.cex = 1.2, which = 1:2)
addtext(g6, c(-0.9, -0.5), c(-0.2, -0.1), "Data Granulo", plabels.cex = 1.2, which = 1:2, plabels.col = c("red", "blue"))

xy2 <- cbind.data.frame(x = runif(200, -1, 1), y = runif(200, -2, 2), y2 = runif(200, -0.5, 0.5))
fac2 <- factor(xy2$x > 0) 
g9 <- s.class(xy2, fac = fac2, xax = 1, yax = 2:3, plot = FALSE)
addtext(g9, 0, 0, "A", plabels.col = "red", plabels.cex = 2)
addtext(g9, c(-2.1, -1.07), c(2, 1), c("A", "B"), plabels.col ="red", plabels.cex = 2, which = 1:2)
addtext(g9, c(-2.1, -1.07), c(2, 1), c("A", "B"), plabels.col = c("green", "red"), plabels.cex = c(3, 2), which = 1:2)

xy <- cbind.data.frame(x = runif(200, -1, 1), y = runif(200, -1, 1))
posi <- factor(xy$x > 0) : factor(xy$y > 0)
g10 <- s.class(xy, fac = posi, facets = posi, pellipses.col = 1:4, plabels.cex = 0, plegend.drawKey = FALSE, psub.cex = 0, plot = FALSE)
addtext(g10, c(0.5, 0.5, -0.5, -0.5), c(0.5, -0.5), levels(posi), plabels.cex = 2, plabels.col = 1:4)

# when superposition
# addtext(g3, -0.9, -0.2, "Data Granulo", plabels.cex = 1.2, which = 1) # message d'erreur : unable to find an inherited method for function ‘printSuperpose’ for signature ‘"S2.label", "ADEgS"’
# addtext(g3, -0.5, -0.1, "Data Granulo", plabels.cex = 1.2, which = 2)
# addtext(g3, -0.9, -0.2, "Data Granulo", plabels.cex = 1.2, which = 1:2)
