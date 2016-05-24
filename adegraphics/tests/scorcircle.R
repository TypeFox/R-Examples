library(adegraphics)
pdf("scorcircle.pdf")
data(olympic, package = "ade4")

dudi1 <- ade4::dudi.pca(olympic$tab, scan = FALSE) # a normed PCA
g1 <- s.corcircle(dudi1$co, lab = names(olympic$tab))

g2 <- s.corcircle(dudi1$co, lab = names(olympic$tab), fullcircle = T)
g3 <- s.corcircle(dudi1$co, lab = names(olympic$tab), fullcircle = FALSE)
g4 <- s.corcircle(dudi1$co, lab = names(olympic$tab), pback.col = "red", pbackground.box = FALSE)
