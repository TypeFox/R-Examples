### R code from vignette source 'sids.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: sids.Rnw:31-36
###################################################
owidth <- getOption("width")
options("width"=70)
ow <- getOption("warn")
options("warn"=-1)
.PngNo <- 0


###################################################
### code chunk number 2: afig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
## opar <- par(mar=c(3,3,1,1)+0.1)


###################################################
### code chunk number 3: afigl (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")


###################################################
### code chunk number 4: bfigl (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=file, width = 6.5, height = 5, pointsize = 12, bg = "white")


###################################################
### code chunk number 5: zfig (eval = FALSE)
###################################################
## par(opar)
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 6: zfigl (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 7: sids.Rnw:77-78
###################################################
library(spdep)


###################################################
### code chunk number 8: sids.Rnw:95-101
###################################################
library(sp)
library(maptools)
library(spdep)
nc_file <- system.file("etc/shapes/sids.shp", package="spdep")[1]
llCRS <- CRS("+proj=longlat +datum=NAD27")
nc <- readShapeSpatial(nc_file, ID="FIPSNO", proj4string=llCRS)


###################################################
### code chunk number 9: plotNC1 (eval = FALSE)
###################################################
## plot(nc, axes=TRUE)
## text(coordinates(nc), label=nc$FIPSNO, cex=0.5)


###################################################
### code chunk number 10: sids.Rnw:123-126
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
opar <- par(mar=c(3,3,1,1)+0.1)
plot(nc, axes=TRUE)
text(coordinates(nc), label=nc$FIPSNO, cex=0.5)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 11: sids.Rnw:139-141
###################################################
names(nc)
summary(nc)


###################################################
### code chunk number 12: sids.Rnw:151-157
###################################################
gal_file <- system.file("etc/weights/ncCR85.gal", package="spdep")[1]
ncCR85 <- read.gal(gal_file, region.id=nc$FIPSNO)
ncCR85
gal_file <- system.file("etc/weights/ncCC89.gal", package="spdep")[1]
ncCC89 <- read.gal(gal_file, region.id=nc$FIPSNO)
ncCC89


###################################################
### code chunk number 13: plot-CC89.nb (eval = FALSE)
###################################################
## plot(nc, border="grey")
## plot(ncCC89, coordinates(nc), add=TRUE, col="blue")


###################################################
### code chunk number 14: sids.Rnw:170-173
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
opar <- par(mar=c(3,3,1,1)+0.1)
plot(nc, border="grey")
plot(ncCC89, coordinates(nc), add=TRUE, col="blue")
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 15: sids.Rnw:193-197
###################################################
ncCC89
r.id <- attr(ncCC89, "region.id")
ncCC89[[match("37001", r.id)]]
r.id[ncCC89[[match("37001", r.id)]]]


###################################################
### code chunk number 16: sids.Rnw:207-208
###################################################
as.character(nc$NAME)[card(ncCC89) == 0]


###################################################
### code chunk number 17: sids.Rnw:243-249
###################################################
ch <- choynowski(nc$SID74, nc$BIR74)
nc$ch_pmap_low <- ifelse(ch$type, ch$pmap, NA)
nc$ch_pmap_high <- ifelse(!ch$type, ch$pmap, NA)
prbs <- c(0,.001,.01,.05,.1,1)
nc$high = cut(nc$ch_pmap_high, prbs)
nc$low = cut(nc$ch_pmap_low,prbs )


###################################################
### code chunk number 18: sids.Rnw:251-252 (eval = FALSE)
###################################################
## spplot(nc, c("low", "high"), col.regions=grey.colors(5))


###################################################
### code chunk number 19: sids.Rnw:258-261
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
print(spplot(nc, c("low", "high"), col.regions=grey.colors(5)))
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 20: sids.Rnw:285-289
###################################################
pmap <- probmap(nc$SID74, nc$BIR74)
nc$pmap <- pmap$pmap
brks <- c(0,0.001,0.01,0.025,0.05,0.95,0.975,0.99,0.999,1)
library(RColorBrewer)


###################################################
### code chunk number 21: sids.Rnw:291-292 (eval = FALSE)
###################################################
## spplot(nc, "pmap", at=brks, col.regions=rev(brewer.pal(9, "RdBu")))


###################################################
### code chunk number 22: sids.Rnw:298-301
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
print(spplot(nc, "pmap", at=brks, col.regions=rev(brewer.pal(9, "RdBu"))))
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 23: poishist (eval = FALSE)
###################################################
## hist(nc$pmap, main="")


###################################################
### code chunk number 24: sids.Rnw:328-331
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
opar <- par(mar=c(3,3,1,1)+0.1)
hist(nc$pmap, main="")
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 25: sids.Rnw:346-349
###################################################
res <- glm(SID74 ~ offset(log(BIR74)), data=nc, family="quasipoisson")
nc$stdres <- rstandard(res)
brks <- c(-4, -3, -2, -1.5, -1, -0.5, 0.5, 1, 1.5, 2, 3, 4)


###################################################
### code chunk number 26: sids.Rnw:351-352 (eval = FALSE)
###################################################
## spplot(nc, "stdres", at=brks, col.regions=rev(brewer.pal(11, "RdBu")))


###################################################
### code chunk number 27: sids.Rnw:370-373
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
print(spplot(nc, "stdres", at=brks, col.regions=rev(brewer.pal(11, "RdBu"))))
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 28: sids.Rnw:402-408
###################################################
global_rate <- sum(nc$SID74)/sum(nc$BIR74)
nc$Expected <- global_rate * nc$BIR74
res <- EBlocal(nc$SID74, nc$Expected, ncCC89, zero.policy=TRUE)
nc$EB_loc <- res$est
brks <- c(0, 0.25, 0.5, 0.75, 1, 2, 3, 4, 5)
spl <- list("sp.text", loc=coordinates(nc)[card(ncCC89) == 0,], txt=rep("*", 2), cex=1.2)


###################################################
### code chunk number 29: sids.Rnw:410-411 (eval = FALSE)
###################################################
## spplot(nc, "EB_loc", at=brks, col.regions=rev(brewer.pal(8, "RdBu")), sp.layout=spl)


###################################################
### code chunk number 30: sids.Rnw:423-426
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
print(spplot(nc, "EB_loc", at=brks, col.regions=rev(brewer.pal(8, "RdBu")), sp.layout=spl))
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 31: sids.Rnw:438-440
###################################################
set.seed(1)
EBImoran.mc(nc$SID74, nc$BIR74, nb2listw(ncCC89, style="B", zero.policy=TRUE), nsim=999, zero.policy=TRUE)


###################################################
### code chunk number 32: sids.Rnw:454-456
###################################################
nc$both <- factor(paste(nc$L_id, nc$M_id, sep=":"))
nboth <- length(table(unclass(nc$both)))


###################################################
### code chunk number 33: sids.Rnw:458-459 (eval = FALSE)
###################################################
## spplot(nc, "both", col.regions=sample(rainbow(nboth)))


###################################################
### code chunk number 34: sids.Rnw:465-468
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
print(spplot(nc, "both", col.regions=sample(rainbow(nboth))))
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 35: sids.Rnw:478-480
###################################################
nc$ft.SID74 <- sqrt(1000)*(sqrt(nc$SID74/nc$BIR74) + sqrt((nc$SID74+1)/nc$BIR74))
stem(round(nc$ft.SID74, 1), scale=2)


###################################################
### code chunk number 36: sids.Rnw:497-499
###################################################
mBIR74 <- tapply(nc$BIR74, nc$both, sum)
mSID74 <- tapply(nc$SID74, nc$both, sum)


###################################################
### code chunk number 37: sids.Rnw:511-515
###################################################
mFT <- sqrt(1000)*(sqrt(mSID74/mBIR74) + sqrt((mSID74+1)/mBIR74))
mFT1 <- t(matrix(mFT, 4, 4, byrow=TRUE))
med <- medpolish(mFT1, na.rm=TRUE, trace.iter=FALSE)
med


###################################################
### code chunk number 38: sids.Rnw:527-531
###################################################
mL_id <- model.matrix(~ as.factor(nc$L_id) -1)
mM_id <- model.matrix(~ as.factor(nc$M_id) -1)
nc$pred <- c(med$overall + mL_id %*% med$row + mM_id %*% med$col)
nc$mp_resid <- nc$ft.SID74 - nc$pred


###################################################
### code chunk number 39: medpolfig (eval = FALSE)
###################################################
## cI_ft <- pretty(nc$ft.SID74, n=9)
## pal_ft <- colorRampPalette(brewer.pal(6, "YlOrBr"))(length(cI_ft)-1)
## p1 <- spplot(nc, c("ft.SID74"), col.regions=pal_ft, at=cI_ft, col="grey30", main="FT transformed SIDS rate")
## p2 <- spplot(nc, c("pred"), col.regions=pal_ft, at=cI_ft, col="grey30", main="Median-polish fit")
## atn <- pretty(nc$mp_resid[nc$mp_resid < 0])
## atp <- pretty(nc$mp_resid[nc$mp_resid >= 0])
## pal <- c(rev(brewer.pal(length(atn-1), "YlOrRd")), brewer.pal(length(atp[-1]), "YlGnBu")[-1])
## p3 <- spplot(nc, "mp_resid", at=c(atn, atp[-1]), col.regions=pal, col="grey30", main="Median-polish residuals")
## plot(p1, split=c(1,1,1,3), more=TRUE)
## plot(p2, split=c(1,2,1,3), more=TRUE)
## plot(p3, split=c(1,3,1,3), more=FALSE)


###################################################
### code chunk number 40: sids.Rnw:550-553
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 5, pointsize = 12, bg = "white")
cI_ft <- pretty(nc$ft.SID74, n=9)
pal_ft <- colorRampPalette(brewer.pal(6, "YlOrBr"))(length(cI_ft)-1)
p1 <- spplot(nc, c("ft.SID74"), col.regions=pal_ft, at=cI_ft, col="grey30", main="FT transformed SIDS rate")
p2 <- spplot(nc, c("pred"), col.regions=pal_ft, at=cI_ft, col="grey30", main="Median-polish fit")
atn <- pretty(nc$mp_resid[nc$mp_resid < 0])
atp <- pretty(nc$mp_resid[nc$mp_resid >= 0])
pal <- c(rev(brewer.pal(length(atn-1), "YlOrRd")), brewer.pal(length(atp[-1]), "YlGnBu")[-1])
p3 <- spplot(nc, "mp_resid", at=c(atn, atp[-1]), col.regions=pal, col="grey30", main="Median-polish residuals")
plot(p1, split=c(1,1,1,3), more=TRUE)
plot(p2, split=c(1,2,1,3), more=TRUE)
plot(p3, split=c(1,3,1,3), more=FALSE)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 41: sids.Rnw:579-587
###################################################
sids.nhbr30.dist <- nbdists(ncCC89, cbind(nc$east, nc$north))
sids.nhbr <- listw2sn(nb2listw(ncCC89, glist=sids.nhbr30.dist, style="B", zero.policy=TRUE))
dij <- sids.nhbr[,3]
n <- nc$BIR74
el1 <- min(dij)/dij
el2 <- sqrt(n[sids.nhbr$to]/n[sids.nhbr$from])
sids.nhbr$weights <- el1*el2
sids.nhbr.listw <- sn2listw(sids.nhbr)


###################################################
### code chunk number 42: sids.Rnw:594-595
###################################################
nc$ft.NWBIR74 <- sqrt(1000)*(sqrt(nc$NWBIR74/nc$BIR74) + sqrt((nc$NWBIR74+1)/nc$BIR74))


###################################################
### code chunk number 43: sids.Rnw:602-609
###################################################
lm_nc <- lm(ft.SID74 ~ 1, data=nc)
outl <- which.max(rstandard(lm_nc))
as.character(nc$names[outl])
W <- listw2mat(sids.nhbr.listw)
W.4 <- W[-outl, -outl]
sids.nhbr.listw.4 <- mat2listw(W.4)
nc2 <- nc[!(1:length(nc$CNTY_ID) %in% outl),]


###################################################
### code chunk number 44: sids.Rnw:616-618
###################################################
ecarIaw <- spautolm(ft.SID74 ~ 1, data=nc2, listw=sids.nhbr.listw.4, weights=BIR74, family="CAR")
summary(ecarIaw)


###################################################
### code chunk number 45: sids.Rnw:625-627
###################################################
ecarIIaw <- spautolm(ft.SID74 ~ both - 1, data=nc2, listw=sids.nhbr.listw.4, weights=BIR74, family="CAR")
summary(ecarIIaw)


###################################################
### code chunk number 46: sids.Rnw:634-636
###################################################
ecarIVaw <- spautolm(ft.SID74 ~ ft.NWBIR74, data=nc2, listw=sids.nhbr.listw.4, weights=BIR74, family="CAR")
summary(ecarIVaw)


###################################################
### code chunk number 47: sids.Rnw:641-642
###################################################
nc2$fitIV <- fitted.values(ecarIVaw)


###################################################
### code chunk number 48: sids.Rnw:644-645 (eval = FALSE)
###################################################
## spplot(nc2, "fitIV", cuts=12, col.regions=grey.colors(13, 0.9, 0.3))


###################################################
### code chunk number 49: sids.Rnw:651-654
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
print(spplot(nc2, "fitIV", cuts=12, col.regions=grey.colors(13, 0.9, 0.3)))
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 50: sids.Rnw:662-663
###################################################
ecarIawll <- spautolm(ft.SID74 ~ 1, data=nc2, listw=sids.nhbr.listw.4, weights=BIR74, family="CAR", llprof=seq(-0.1, 0.9020532358, length.out=100))


###################################################
### code chunk number 51: ll_prof (eval = FALSE)
###################################################
## plot(ll ~ lambda, ecarIawll$llprof, type="l")


###################################################
### code chunk number 52: sids.Rnw:672-675
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
opar <- par(mar=c(3,3,1,1)+0.1)
plot(ll ~ lambda, ecarIawll$llprof, type="l")
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 53: sids.Rnw:683-685
###################################################
options("width"=owidth)
options("warn"=ow)


