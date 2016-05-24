library(adegraphics)
pdf("parameter.pdf")

adegparold <- adegpar()
b1 <- length(adegraphics:::separation(plines = list(col = "blue"), plab.bo.dra = FALSE, plab = list(orien = F))$rest) == 0
b2 <- length(adegraphics:::separation(plines = list(col = "blue", lwd = c(1:5)), parr.end = NA, plab.boxes.dr = FALSE)$rest) == 0
b3 <- length(adegraphics:::separation(plot.li = list(col = "blue", lwd = c(1:5)), par.end = NA, pl.boxes.draw = FALSE, pattern = 1)$rest) == 2
b4 <- length(adegraphics:::separation(plot.li = list(col = "blue", lwd = c(1:5)), par.end = NA, pl.boxes.draw = FALSE, pattern = 0)$rest) == 2
b5 <- length(adegraphics:::separation()) == 2
b6 <- length(adegraphics:::separation(par.sub.text = list("lineheight" = 5), pattern = 1)$rest) == 0
b7 <- names(adegraphics:::separation(toto = "rienavoir!")$rest) == "toto"
b8 <- names(adegraphics:::separation(toto = "rienavoir!", pattern = 1)$rest) == "toto"

l1 <- list(parrow = list(col = "blue", lwd = c(1:5), end = NA), plboxes.draw = FALSE, pattern = 1)
sep1 <- adegraphics:::separation(l1)  ## no recognition of "pattern" in a list

l2 <- list("parrow.lwd", "plboxes.draw")
sep2 <- adegraphics:::separation(l2)

sep3 <- adegraphics:::separation("plboxes" = list("draw" = FALSE, "col"))
sep4 <- adegraphics:::separation("plboxes" = list("draw" = FALSE, "col" = 2))
sep5 <- adegraphics:::separation(pla.box = list(col = 1:5))
sep6 <- adegraphics:::separation(pla.box.col = c(1:5))
sep7 <- adegraphics:::separation(pla = list(box.col = 1:5))  ## don't match


## adegpar test
ad1 <- adegpar()
ad2 <- adegpar("paxes.draw")
ad3 <- adegpar("paxes.draw", "psub.cex")

ad4 <- adegpar("psub.cex" = 5)
ad5 <- adegpar("psub.cex" = 5, paxes.draw = FALSE)
ad6 <- adegpar("psub")
ad7 <- adegpar("psub.cex", "plabe.boxes")

ad8 <- adegpar(ppoints = list(col = "yellow"), pgrid.space = 4, plines = list(lwd = c(1:5)))
ad9 <- adegpar(ppoints = list(col = "red"), pgrid = list(nint = 12), plines = list(lwd = c(1:5)))
ad10 <- adegpar("ppoints.col", "pgrid.nint", "plines")

ad11 <- adegpar(paxes = list("x"), pgrid = list("nint", "col"))
ad12 <- adegpar(list(pellip = list(col = "red"), grid.nint = 8, plines = list(lwd = c(1:5))))
ad13 <- adegpar(paxes = list("x"), pgrid = list("nint", "col"), "plines", "pellipse")

ad14 <- adegpar(plegend.drawKey = FALSE)
ad15 <- adegpar(list(paxes = list(col = "white"), pgrid.nint = 6, plines = list(lwd = c(1:5))))
ad16 <- adegpar(paxes = list(x = list(draw = TRUE)))

adegpar(adegparold)

## merging list
l3 <- list(plabels = list(boxes = list(col = "white", alpha = 1)), plabels = list(cex = 2), plabels = list(col = "red"))
adegraphics:::.mergingList(l3)
adegraphics:::.mergingList(list(plabels = list(cex = 3)))


## update parameters in graphics
cha <- rep(LETTERS, length.out = 100)
xy <- cbind.data.frame(runif(length(cha)), runif(length(cha)))
g1 <- s.label(xy, labels = cha, paxes.draw = TRUE, plabels.cex = runif(length(cha), 0.5, 1.5))
update(g1, paxes = list(aspect = "fill", draw = TRUE, x = list(draw = FALSE)), 
  pgrid = list(col = "black", lwd = 2, lty = 5), 
  plabels = list(col = 1:4, alpha = 0.5, cex = 2, boxes = list(border = "blue", col = "antiquewhite", alpha = 0.2, lwd = 2.5, lty = 5))) 

g2 <- s.label(xy, labels = cha)
update(g2, paxes.draw = FALSE, pgrid = list(col = "blue", lwd = 2, lty = 5, text = list(pos = "bottomright", cex = 2, col = "green")))
update(g2, porigin = list(alpha = 0.5, col = "red", lty = 5 , lwd = 2, origin = c(0.6, 0.1)), pgrid.lwd = 0.5)
update(g2, psub = list(text = "parameters", cex = 2, col = "red", position = "topright"))
update(g2, plabels.cex = 0, ppoints = list(alpha = 0.4, cex = 2, col = "red", fill = "blue", pch = 21))


## from tdr641
data(doubs, package = "ade4")
dudi1 <- ade4::dudi.pca(doubs$env, scale = T, scannf = F, nf = 3)
dudi2 <- ade4::dudi.pca(doubs$fish, scale = T, scannf = F, nf = 2)
coin1 <- ade4::coinertia(dudi1, dudi2, scannf = F, nf = 2)
g3 <- s.arrow(coin1$l1, plabels.cex = .87)
update(g3, plines = list(col = "blue", lwd = 2, lty = 3), parr.end = "both", parr = list(angle = 25, length = 0.5))

## with spatial object
library(maptools)
nc <- readShapePoly(system.file("shapes/sids.shp", package = "maptools")[1], proj4string = CRS("+proj=longlat +datum=NAD27"))
xy <- coordinates(nc)
g4 <- s.label(xy, label = as.character(1:nrow(xy)), porigin.include = FALSE, Sp = nc, pSp.col = colorRampPalette(c("yellow", "blue"))(52), pgrid.draw = TRUE)
update(g4, pSp = list(col = "yellow", border = "blue", lwd = 2, lty = 5, alpha = 0.01)) ## don't match : to solve

## plabels parameter
data(tortues, package = "ade4")
pturtles <- tortues
names(pturtles) <- c("length", "width", "height", "sex")
sex <- pturtles$sex
sexcol <- ifelse(sex == "M", "blue", "red")
measures <- pturtles[, 1:3]
pca1 <- ade4::dudi.pca(measures, scann = FALSE, nf = 3)
g5 <- scatter(pca1, row.plabel.cex = 0, col.plabel.cex = c(1, 2, 3), posieig = "none", col.plabel.col = c("red", "blue", "green"))
