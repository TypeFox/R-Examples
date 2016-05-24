library(adegraphics)
pdf("shist.pdf")

dfxy1 <- matrix(rnorm(200), ncol = 2)
g1 <- s.label(dfxy1)
g2 <- addhist(g1)

dfxy2 <- dfxy1
dfxy2[, 2] <- dfxy2[, 2] + rnorm(100, 2)
g3 <- s.label(dfxy2)
g4 <- addhist(g3, plot.polygon = list(col = "red"))

data(rpjdl, package = "ade4")
coa1 <- ade4::dudi.coa(rpjdl$fau, scannf = FALSE, nf = 4)
g5 <- s.label(coa1$li)
g6 <- addhist(g5)
