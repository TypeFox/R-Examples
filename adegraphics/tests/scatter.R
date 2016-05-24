## delete/remove this file when 'scatter' functions will be removed in ade4

library(adegraphics)
pdf("scatter.pdf")

data(deug, package = "ade4")
dd1 <- ade4::dudi.pca(deug$tab, scannf = FALSE, nf = 4)
scatter(dd1, posieig = "bottomright")
scatter(dd1, posieig = "bottomright", plot = T, prop = TRUE)
scatter(dd1, posieig = "none", plot = T)
scatter(dd1, posieig = "bottomleft", plot = T)
scatter(dd1, posieig = "topright", plot = T)
scatter(dd1, posieig = "topleft", plot = T, eig.col = c("white", "blue", "red"))

data(rhone, package = "ade4")
dd1 <- ade4::dudi.pca(rhone$tab, nf = 4, scannf = FALSE)
g1 <- scatter(dd1, sub = "Principal component analysis", row = list(plabels.optim = TRUE), col.pla.boxes.alpha = 0.5)
g1[2, drop = TRUE]
scatter(dd1, row = list(sub = "Principal component analysis", plabels.optim = TRUE), col.pla.boxes.alpha = 0.5)
scatter(dd1, prop = TRUE, ppoints.cex = 0.2, density.plot = TRUE, row = list(threshold = 0.01))


##################### scatter.coa test
data(housetasks, package = "ade4")
par(mfrow = c(2, 2))
dd2 <- ade4::dudi.coa(housetasks, scan = FALSE)
ade4::scatter(dd2, method = 1, sub = "1 / Standard", posieig = "none")
ade4::scatter(dd2, method = 2, sub = "2 / Columns -> averaging -> Rows", posieig = "none")
ade4::scatter(dd2, method = 3, sub = "3 / Rows -> averaging -> Columns ", posieig = "none")
g1 <- scatter(dd2, method = 1, row.sub = "1 / Standard", posieig = "none", plot = FALSE)
g2 <- scatter(dd2, method = 2, col.sub = "2 / Columns -> averaging -> Rows", posieig = "none", plot = FALSE)
g3 <- scatter(dd2, method = 3, row.sub = "3 / Rows -> averaging -> Columns ", posieig = "none", plot = FALSE)
G <- ADEgS(list(g1, g2, g3), layout = c(2, 2), plot = TRUE)


##################### plot.acm test
data(lascaux, package = "ade4")

acm1 <- ade4::dudi.acm(lascaux$ornem, sca = FALSE)
p1 <- proc.time()
ade4::scatter(acm1)
Tade4 <- proc.time() - p1

p2 <- proc.time()
plot(acm1, ppoints.cex = 0.3, plot = T)
Tadegraphics <- proc.time() - p2
## faster caculus, longest display than for ade4


##################### plot.fca text
data(coleo, package = "ade4")
coleo.fuzzy <- ade4::prep.fuzzy.var(coleo$tab, coleo$col.blocks)

fca1 <- ade4::dudi.fca(coleo.fuzzy, scannf = FALSE, nf = 3)
ade4::scatter(fca1)
plot(fca1)
