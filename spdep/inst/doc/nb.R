### R code from vignette source 'nb.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: nb.Rnw:43-48
###################################################
owidth <- getOption("width")
options("width"=90)
ow <- getOption("warn")
options("warn"=-1)
.PngNo <- 0


###################################################
### code chunk number 2: figreset (eval = FALSE)
###################################################
## .iwidth <- 5
## .iheight <- 6
## .ipointsize <- 12


###################################################
### code chunk number 3: nb.Rnw:56-57
###################################################
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 4: afig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
## opar <- par(mar=c(3,3,1,1)+0.1)


###################################################
### code chunk number 5: afigl (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")


###################################################
### code chunk number 6: bfigl (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=file, width = 6.5, height = 5, pointsize = 12, bg = "white")


###################################################
### code chunk number 7: bfig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=file, width = 6.5, height = 5, pointsize = 12, bg = "white")
## opar <- par(mar=c(3,3,1,1)+0.1)


###################################################
### code chunk number 8: zfig (eval = FALSE)
###################################################
## par(opar)
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 9: zfigl (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 10: nb.Rnw:113-121
###################################################
library(maptools)
owd <- getwd()
setwd(system.file("etc/shapes", package="spdep"))
NY8 <- readShapeSpatial("NY8_utm18")
library(spdep)
setwd(system.file("etc/weights", package="spdep"))
NY_nb <- read.gal("NY_nb.gal", region.id=row.names(NY8))
setwd(owd)


###################################################
### code chunk number 11: nb.Rnw:133-136
###################################################
Syracuse <- NY8[NY8$AREANAME == "Syracuse city",]
Sy0_nb <- subset(NY_nb, NY8$AREANAME == "Syracuse city")
summary(Sy0_nb)


###################################################
### code chunk number 12: nb.Rnw:152-155
###################################################
class(Syracuse)
Sy1_nb <- poly2nb(Syracuse)
isTRUE(all.equal(Sy0_nb, Sy1_nb, check.attributes=FALSE))


###################################################
### code chunk number 13: nb.Rnw:180-182
###################################################
Sy2_nb <- poly2nb(Syracuse, queen=FALSE)
isTRUE(all.equal(Sy0_nb, Sy2_nb, check.attributes=FALSE))


###################################################
### code chunk number 14: nb.Rnw:189-205
###################################################
.iwidth <- 6
.iheight <- 4.5
.ipointsize <- 10
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
opar <- par(mar=c(3,3,1,1)+0.1)
oopar <- par(mfrow=c(1,2), mar=c(3,3,1,1)+0.1)
plot(Syracuse, border="grey60")
plot(Sy0_nb, coordinates(Syracuse), add=TRUE, pch=19, cex=0.6)
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="a)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy0_nb, coordinates(Syracuse), add=TRUE, pch=19, cex=0.6)
plot(diffnb(Sy0_nb, Sy2_nb, verbose=FALSE), coordinates(Syracuse),
  add=TRUE, pch=".", cex=0.6, lwd=2)
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="b)", cex=0.8)
par(oopar)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 15: nb.Rnw:239-242 (eval = FALSE)
###################################################
## library(spgrass6)
## writeVECT6(Syracuse, "SY0")
## contig <- vect2neigh("SY0")


###################################################
### code chunk number 16: nb.Rnw:244-245 (eval = FALSE)
###################################################
## load("contig.RData")


###################################################
### code chunk number 17: nb.Rnw:247-249 (eval = FALSE)
###################################################
## Sy3_nb <- sn2listw(contig)$neighbours
## isTRUE(all.equal(Sy3_nb, Sy2_nb, check.attributes=FALSE))


###################################################
### code chunk number 18: nb.Rnw:300-307
###################################################
coords <- coordinates(Syracuse)
IDs <- row.names(as(Syracuse, "data.frame"))
#FIXME library(tripack)
Sy4_nb <- tri2nb(coords, row.names=IDs)
Sy5_nb <- graph2nb(soi.graph(Sy4_nb, coords), row.names=IDs)
Sy6_nb <- graph2nb(gabrielneigh(coords), row.names=IDs)
Sy7_nb <- graph2nb(relativeneigh(coords), row.names=IDs)


###################################################
### code chunk number 19: nb.Rnw:313-333
###################################################
.iwidth <- 5
.iheight <- 3.5
.ipointsize <- 10
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
opar <- par(mar=c(3,3,1,1)+0.1)
oopar <- par(mfrow=c(2,2), mar=c(1,1,1,1)+0.1)
plot(Syracuse, border="grey60")
plot(Sy4_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="a)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy5_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="b)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy6_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="c)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy7_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="d)", cex=0.8)
par(oopar)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 20: nb.Rnw:363-367
###################################################
nb_l <- list(Triangulation=Sy4_nb, SOI=Sy5_nb, Gabriel=Sy6_nb,
  Relative=Sy7_nb)
sapply(nb_l, function(x) is.symmetric.nb(x, verbose=FALSE, force=TRUE))
sapply(nb_l, function(x) n.comp.nb(x)$nc)


###################################################
### code chunk number 21: nb.Rnw:386-392
###################################################
Sy8_nb <- knn2nb(knearneigh(coords, k=1), row.names=IDs)
Sy9_nb <- knn2nb(knearneigh(coords, k=2), row.names=IDs)
Sy10_nb <- knn2nb(knearneigh(coords, k=4), row.names=IDs)
nb_l <- list(k1=Sy8_nb, k2=Sy9_nb, k4=Sy10_nb)
sapply(nb_l, function(x) is.symmetric.nb(x, verbose=FALSE, force=TRUE))
sapply(nb_l, function(x) n.comp.nb(x)$nc)


###################################################
### code chunk number 22: nb.Rnw:398-415
###################################################
.iwidth <- 5
.iheight <- 2.5
.ipointsize <- 10
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
opar <- par(mar=c(3,3,1,1)+0.1)
oopar <- par(mfrow=c(1,3), mar=c(1,1,1,1)+0.1)
plot(Syracuse, border="grey60")
plot(Sy8_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="a)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy9_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="b)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy10_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="c)", cex=0.8)
par(oopar)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 23: nb.Rnw:439-449
###################################################
dsts <- unlist(nbdists(Sy8_nb, coords))
summary(dsts)
max_1nn <- max(dsts)
max_1nn
Sy11_nb <- dnearneigh(coords, d1=0, d2=0.75*max_1nn, row.names=IDs)
Sy12_nb <- dnearneigh(coords, d1=0, d2=1*max_1nn, row.names=IDs)
Sy13_nb <- dnearneigh(coords, d1=0, d2=1.5*max_1nn, row.names=IDs)
nb_l <- list(d1=Sy11_nb, d2=Sy12_nb, d3=Sy13_nb)
sapply(nb_l, function(x) is.symmetric.nb(x, verbose=FALSE, force=TRUE))
sapply(nb_l, function(x) n.comp.nb(x)$nc)


###################################################
### code chunk number 24: nb.Rnw:455-472
###################################################
.iwidth <- 5
.iheight <- 2.5
.ipointsize <- 10
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
opar <- par(mar=c(3,3,1,1)+0.1)
oopar <- par(mfrow=c(1,3), mar=c(1,1,1,1)+0.1)
plot(Syracuse, border="grey60")
plot(Sy11_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="a)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy12_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="b)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy13_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="c)", cex=0.8)
par(oopar)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 25: nb.Rnw:481-498
###################################################
.iwidth <- 5
.iheight <- 2.5
.ipointsize <- 10
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
opar <- par(mar=c(3,3,1,1)+0.1)
oopar <- par(mfrow=c(1,3), mar=c(1,1,1,1)+0.1)
plot(Syracuse, border="grey60")
plot(Sy11_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="a)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy12_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="b)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy13_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="c)", cex=0.8)
par(oopar)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")
.iwidth <- 5
.iheight <- 6
.ipointsize <- 12


###################################################
### code chunk number 26: nb.Rnw:523-525
###################################################
dsts0 <- unlist(nbdists(NY_nb, coordinates(NY8)))
summary(dsts0)


###################################################
### code chunk number 27: nb.Rnw:557-558
###################################################
Sy0_nb_lags <- nblag(Sy0_nb, maxlag=9)


###################################################
### code chunk number 28: nb.Rnw:571-581
###################################################
library(xtable)
names(Sy0_nb_lags) <- c("first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth")
res <- sapply(Sy0_nb_lags, function(x) table(card(x)))
mx <- max(unlist(sapply(Sy0_nb_lags, function(x) card(x))))
nn <- length(Sy0_nb_lags)
res1 <- matrix(0, ncol=(mx+1), nrow=nn)
rownames(res1) <- names(res)
colnames(res1) <- as.character(0:mx)
for (i in 1:nn) res1[i, names(res[[i]])] <- res[[i]]
print(xtable(t(res1), align=c("r", "|", rep("r", nn)), digits=0), floating=FALSE)


###################################################
### code chunk number 29: nb.Rnw:608-610
###################################################
cell2nb(7, 7, type="rook", torus=TRUE)
cell2nb(7, 7, type="rook", torus=FALSE)


###################################################
### code chunk number 30: nb.Rnw:627-634
###################################################
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
gridded(meuse.grid) <- TRUE
dst <- max(slot(slot(meuse.grid, "grid"), "cellsize"))
mg_nb <- dnearneigh(coordinates(meuse.grid), 0, dst)
mg_nb
table(card(mg_nb))


