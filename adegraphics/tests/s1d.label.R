library(adegraphics)
pdf("s1d.label.pdf")
data(meau, package= "ade4")

envpca <- ade4::dudi.pca(meau$env, scannf = FALSE)

g1 <- s1d.label(envpca$l1[, 1], row.names(envpca$l1))
g2 <- s1d.label(envpca$l1[, 1], row.names(envpca$l1), p1d.hori = F)
g3 <- s1d.label(envpca$l1[, 1], row.names(envpca$l1), plabels.boxes.draw = FALSE, plab.srt = 45, plabel.boxes =  list(draw = FALSE))
g4 <- s1d.label(envpca$co[, 1], row.names(envpca$co), p1d.reverse = TRUE, poslabel = "value")
g5 <- s1d.label(envpca$l1[, 1], row.names(envpca$l1), at = 0, plabel.cex = 0)
