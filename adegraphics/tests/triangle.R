library(adegraphics)
pdf("triangle.pdf")

## ex1
data(euro123, package = "ade4")

dfxyz1 <- rbind.data.frame(euro123$in78, euro123$in86, euro123$in97)
row.names(dfxyz1) <- paste(row.names(euro123$in78), rep(c(1, 2, 3), rep(12, 3)), sep = "")
g1 <- triangle.label(dfxyz1, label = row.names(dfxyz1))
g2 <- triangle.label(euro123$in86, label = row.names(euro123$in78), plab.cex = 0.8)
g3 <- triangle.match(euro123$in78, euro123$in86)
g4 <- triangle.label(rbind.data.frame(euro123$in78, euro123$in86), plab.cex = 1, addaxes = TRUE, psub = list(text = "Principal axis", cex = 2, pos = "topright"))
g5 <- triangle.label(euro123[[1]], min3 = c(0, 0.2, 0.3), max3 = c(0.5, 0.7, 0.8), plabels.cex = 1, label = row.names(euro123[[1]]), addax = TRUE)
g6 <- triangle.label(euro123[[2]], min3 = c(0, 0.2, 0.3), max3 = c(0.5, 0.7, 0.8), label = row.names(euro123[[1]]), addax = TRUE)
g7 <- triangle.label(euro123[[3]], min3 = c(0, 0.2, 0.3), max3 = c(0.5, 0.7, 0.8), label = row.names(euro123[[1]]), addax = TRUE)
g8 <- triangle.label(rbind.data.frame(euro123[[1]], euro123[[2]], euro123[[3]]))

## ex2
dfxyz2 <- cbind.data.frame(a = runif(100), b = runif(100), c = runif(100, 4, 5))
g9 <- triangle.label(dfxyz2)

## ex3
g10 <- triangle.label(dfxyz1)
g11 <- triangle.class(dfxyz1, as.factor(rep("G", 36)), star = 0.5, ellips = 1)
g12 <- triangle.class(dfxyz1, euro123$plan$an)
g13 <- triangle.class(dfxyz1, euro123$plan$pays)
g14 <- triangle.class(dfxyz1, euro123$plan$an, elli = 1, pell.axe.draw = TRUE)
g15 <- triangle.class(dfxyz1, euro123$plan$an, elli = 0, sta = 0, col = c("red", "green", "blue"), pell.axe.draw = TRUE, plab.cex = 2, ppoi.cex = 2, pell.axe.draw = TRUE)
g16 <- triangle.class(dfxyz1, euro123$plan$an, ell = 2, sta = 0.5, pell.axe.draw = TRUE, plab.cex = 1.5)
g17 <- triangle.class(dfxyz1, euro123$plan$an, ell = 0, sta = 1, adjust = FALSE) 
g18 <- triangle.class(dfxyz1, euro123$plan$an, ell = 0, sta = 1, chull =c(0.2, 0.25, 0.5, 0.75, 1), adjust = TRUE, showposi = TRUE, col = 10:13, pgrid.draw = FALSE)
