### R code from vignette source 'CO69.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: CO69.Rnw:38-43
###################################################
owidth <- getOption("width")
options("width"=90)
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
### code chunk number 5: bfig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
## pdf(file=file, width = 6.5, height = 5, pointsize = 12, bg = "white")
## opar <- par(mar=c(3,3,1,1)+0.1)


###################################################
### code chunk number 6: zfig (eval = FALSE)
###################################################
## par(opar)
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 7: zfigl (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 8: CO69.Rnw:100-105
###################################################
library(spdep)
require(maptools)
fn <- system.file("etc/shapes/eire.shp", package="spdep")[1]
prj <- CRS("+proj=utm +zone=30 +units=km +ellps=mod_airy")
eire <- readShapeSpatial(fn, ID="names", proj4string=prj)


###################################################
### code chunk number 9: CO69.Rnw:107-109
###################################################
class(eire)
names(eire)


###################################################
### code chunk number 10: CO69.Rnw:122-125
###################################################
fn <- system.file("etc/misc/geary_eire.txt", package="spdep")[1]
ge <- read.table(fn, header=TRUE)
names(ge)


###################################################
### code chunk number 11: CO69.Rnw:133-136
###################################################
row.names(ge) <- as.character(ge$county)
all.equal(row.names(ge), row.names(eire))
eire_ge <- spCbind(eire, ge)


###################################################
### code chunk number 12: CO69.Rnw:144-146
###################################################
eire_ge1 <- eire_ge[!(row.names(eire_ge) %in% "Dublin"),]
length(row.names(eire_ge1))


###################################################
### code chunk number 13: CO69.Rnw:153-155
###################################################
skewness <- function(z) {z <- scale(z, scale=FALSE); ((sum(z^3)/length(z))^2)/((sum(z^2)/length(z))^3)}
kurtosis <- function(z) {z <- scale(z, scale=FALSE); (sum(z^4)/length(z))/((sum(z^2)/length(z))^2)}


###################################################
### code chunk number 14: CO69.Rnw:162-166
###################################################
print(sapply(as(eire_ge1, "data.frame")[13:24], skewness), digits=3)
print(sapply(as(eire_ge1, "data.frame")[13:24], kurtosis), digits=4)
print(sapply(as(eire_ge1, "data.frame")[c(13,16,18,19)], function(x) skewness(log(x))), digits=3)
print(sapply(as(eire_ge1, "data.frame")[c(13,16,18,19)], function(x) kurtosis(log(x))), digits=4)


###################################################
### code chunk number 15: CO69.Rnw:212-215
###################################################
fn <- system.file("etc/misc/unstand_sn.txt", package="spdep")[1]
unstand <- read.table(fn, header=TRUE)
summary(unstand)


###################################################
### code chunk number 16: CO69.Rnw:225-231
###################################################
class(unstand) <- c("spatial.neighbour", class(unstand))
of <- ordered(unstand$from)
attr(unstand, "region.id") <- levels(of)
unstand$from <- as.integer(of)
unstand$to <- as.integer(ordered(unstand$to))
attr(unstand, "n") <- length(unique(unstand$from))


###################################################
### code chunk number 17: CO69.Rnw:240-243
###################################################
lw_unstand <- sn2listw(unstand)
lw_unstand$style <- "B"
lw_unstand


###################################################
### code chunk number 18: CO69.Rnw:257-259
###################################################
nb <- poly2nb(eire_ge1)
nb


###################################################
### code chunk number 19: CO69.Rnw:268-269
###################################################
xx <- diffnb(nb, lw_unstand$neighbours, verbose=TRUE)


###################################################
### code chunk number 20: CO69.Rnw:271-274 (eval = FALSE)
###################################################
## plot(eire_ge1, border="grey60")
## plot(nb, coordinates(eire_ge1), add=TRUE, pch=".", lwd=2)
## plot(xx, coordinates(eire_ge1), add=TRUE, pch=".", lwd=2, col=3)


###################################################
### code chunk number 21: CO69.Rnw:280-292
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
opar <- par(mar=c(3,3,1,1)+0.1)
par(mfrow=c(1,2))
plot(eire_ge1, border="grey40")
title(xlab="25 Irish counties")
text(coordinates(eire_ge1), labels=as.character(eire_ge1$serlet), cex=0.8)
plot(eire_ge1, border="grey60")
title(xlab="Contiguities")
plot(nb, coordinates(eire_ge1), add=TRUE, pch=".", lwd=2)
plot(xx, coordinates(eire_ge1), add=TRUE, pch=".", lwd=2, col=3)
legend("topleft", legend=c("Contiguous", "Ferry"), lwd=2, lty=1, col=c(1,3), bty="n", cex=0.7)
par(mfrow=c(1,1))
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 22: CO69.Rnw:322-323
###################################################
load(system.file("etc/misc/raw_grass_borders.RData", package="spdep")[1])


###################################################
### code chunk number 23: CO69.Rnw:325-332 (eval = FALSE)
###################################################
## library(maptools)
## SG <- Sobj_SpatialGrid(eire_ge1)$SG
## library(spgrass6)
## grass_home <- "/home/rsb/topics/grass/g64/grass-6.4.0svn"
## initGRASS(grass_home, home=tempdir(), SG=SG, override=TRUE)
## writeVECT6(eire_ge1, "eire", v.in.ogr_flags=c("o", "overwrite"))
## res <- vect2neigh("eire", ID="serlet")


###################################################
### code chunk number 24: CO69.Rnw:334-343
###################################################
grass_borders <- sn2listw(res)
raw_borders <- grass_borders$weights
int_tot <- attr(res, "total") - attr(res, "external")
prop_borders <- lapply(1:length(int_tot), function(i) raw_borders[[i]]/int_tot[i])
dlist <- nbdists(grass_borders$neighbours, coordinates(eire_ge1))
inv_dlist <- lapply(dlist, function(x) 1/(x/1.609344))
combo_km <- lapply(1:length(inv_dlist), function(i) inv_dlist[[i]]*prop_borders[[i]])
combo_km_lw <- nb2listw(grass_borders$neighbours, glist=combo_km, style="B")
summary(combo_km_lw)


###################################################
### code chunk number 25: CO69.Rnw:356-367
###################################################
red_lw_unstand <- lw_unstand
Clare <- which(attr(lw_unstand, "region.id") == "C")
Kerry <- which(attr(lw_unstand, "region.id") == "H")
Kerry_in_Clare <- which(lw_unstand$neighbours[[Clare]] == Kerry)
Clare_in_Kerry <- which(lw_unstand$neighbours[[Kerry]] == Clare)
red_lw_unstand$neighbours[[Clare]] <- red_lw_unstand$neighbours[[Clare]][-Kerry_in_Clare]
red_lw_unstand$neighbours[[Kerry]] <- red_lw_unstand$neighbours[[Kerry]][-Clare_in_Kerry]
red_lw_unstand$weights[[Clare]] <- red_lw_unstand$weights[[Clare]][-Kerry_in_Clare]
red_lw_unstand$weights[[Kerry]] <- red_lw_unstand$weights[[Kerry]][-Clare_in_Kerry]
summary(red_lw_unstand)
cor(unlist(red_lw_unstand$weights), unlist(combo_km_lw$weights))


###################################################
### code chunk number 26: CO69.Rnw:380-392
###################################################
flatten <- function(x, digits=3, statistic="I") {
  res <- c(format(x$estimate, digits=digits),
    format(x$statistic, digits=digits),
    format.pval(x$p.value, digits=digits))
  res <- matrix(res, ncol=length(res))
  colnames(res) <- paste(c("", "E", "V", "SD_", "P_"), "I", sep="")
  rownames(res) <- deparse(substitute(x))
  res
}
`reconstructed weights` <- moran.test(eire_ge1$ocattlepacre, combo_km_lw)
`original weights` <- moran.test(eire_ge1$ocattlepacre, red_lw_unstand)
print(rbind(flatten(`reconstructed weights`), flatten(`original weights`)), quote=FALSE)


###################################################
### code chunk number 27: CO69.Rnw:413-425
###################################################
eire_ge1$ln_pagval2_10 <- log(eire_ge1$pagval2_10)
eire_ge1$ln_cowspacre <- log(eire_ge1$cowspacre)
eire_ge1$ln_pigspacre <- log(eire_ge1$pigspacre)
eire_ge1$ln_sheeppacre <- log(eire_ge1$sheeppacre)
vars <- c("pagval2_10", "ln_pagval2_10", "pagval10_50", "pagval50p",
 "cowspacre", "ln_cowspacre", "ocattlepacre", "pigspacre",
 "ln_pigspacre", "sheeppacre", "ln_sheeppacre", "townvillp",
 "carspcap", "radiopcap", "retailpcap", "psinglem30_34")
nb_B <- nb2listw(lw_unstand$neighbours, style="B")
nb_B
lw_std <- nb2listw(lw_unstand$neighbours, glist=lw_unstand$weights, style="W")
lw_std


###################################################
### code chunk number 28: CO69.Rnw:455-468
###################################################
system.time({
MoranN <- lapply(vars, function(x) moran.test(eire_ge1[[x]], listw=nb_B, randomisation=FALSE))
MoranR <- lapply(vars, function(x) moran.test(eire_ge1[[x]], listw=nb_B, randomisation=TRUE))
GearyN <- lapply(vars, function(x) geary.test(eire_ge1[[x]], listw=nb_B, randomisation=FALSE))
GearyR <- lapply(vars, function(x) geary.test(eire_ge1[[x]], listw=nb_B, randomisation=TRUE))
Prop_unstdN  <- lapply(vars, function(x) moran.test(eire_ge1[[x]], listw=lw_unstand, randomisation=FALSE))
Prop_unstdR  <- lapply(vars, function(x) moran.test(eire_ge1[[x]], listw=lw_unstand, randomisation=TRUE))
Prop_stdN  <- lapply(vars, function(x) moran.test(eire_ge1[[x]], listw=lw_std, randomisation=FALSE))
Prop_stdR  <- lapply(vars, function(x) moran.test(eire_ge1[[x]], listw=lw_std, randomisation=TRUE))
})
res <- sapply(c("MoranN", "MoranR", "GearyN", "GearyR", "Prop_unstdN", "Prop_unstdR", "Prop_stdN", "Prop_stdR"), function(x) sapply(get(x), "[[", "statistic"))
rownames(res) <- vars
ores <- res[,c(1,2,5:8)]


###################################################
### code chunk number 29: CO69.Rnw:480-481
###################################################
options("width"=100)


###################################################
### code chunk number 30: CO69.Rnw:483-484
###################################################
print(formatC(res, format="f", digits=4), quote=FALSE)


###################################################
### code chunk number 31: CO69.Rnw:486-487
###################################################
options("width"=90)


###################################################
### code chunk number 32: CO69.Rnw:499-507
###################################################
wc_unstd <- spweights.constants(lw_unstand)
wrong_N_sqVI <- sqrt((wc_unstd$nn*wc_unstd$S1 - wc_unstd$n*wc_unstd$S2 + 3*wc_unstd$S0*wc_unstd$S0)/((wc_unstd$nn-1)*wc_unstd$S0*wc_unstd$S0)-((-1/(wc_unstd$n-1))^2))
raw_data <- grep("^ln_", vars, invert=TRUE)
I <- sapply(Prop_stdN, function(x) x$estimate[1])[raw_data]
EI <- sapply(Prop_stdN, function(x) x$estimate[2])[raw_data]
res <- (I - EI)/wrong_N_sqVI
names(res) <- vars[raw_data]
print(formatC(res, format="f", digits=4), quote=FALSE)


###################################################
### code chunk number 33: CO69.Rnw:517-522
###################################################
res <- lapply(c("MoranR", "GearyR", "Prop_unstdR", "Prop_stdR"), function(x) sapply(get(x), function(y) c(y$estimate[1], sqrt(y$estimate[3]))))
res <- t(do.call("rbind", res))
colnames(res) <- c("I", "sigma_I", "C", "sigma_C", "unstd_r", "sigma_r", "std_r", "sigma_r")
rownames(res) <- vars
print(formatC(res, format="f", digits=4), quote=FALSE)


###################################################
### code chunk number 34: CO69.Rnw:539-550
###################################################
oMoranf <- function(x, nb) {
  z <- scale(x, scale=FALSE)
  n <- length(z)
  glist <- lapply(1:n, function(i) {ii <- nb[[i]]; ifelse(ii > i, 1, 0)})
  lw <- nb2listw(nb, glist=glist, style="B")
  wz <- lag(lw, z)
  I <- (sum(z*wz)/sum(z*z))
  I
}
res <- sapply(vars, function(x) oMoranf(eire_ge1[[x]], nb=lw_unstand$neighbours))
print(formatC(res, format="f", digits=4), quote=FALSE)


###################################################
### code chunk number 35: CO69.Rnw:582-601
###################################################
MoranI.boot <- function(var, i, ...) {
  var <- var[i]
  return(moran(x=var, ...)$I)
}
Nsim <- function(d, mle) {
  n <- length(d)
  rnorm(n, mle$mean, mle$sd)
}
f_bperm <- function(x, nsim, listw) {
  boot(x, statistic=MoranI.boot, R=nsim, sim="permutation", listw=listw,
    n=length(x), S0=Szero(listw))
}
f_bpara <- function(x, nsim, listw) {
  boot(x, statistic=MoranI.boot, R=nsim, sim="parametric", ran.gen=Nsim,
    mle=list(mean=mean(x), sd=sd(x)), listw=listw, n=length(x),
    S0=Szero(listw))
}
nsim <- 4999
set.seed(1234)


###################################################
### code chunk number 36: CO69.Rnw:616-624 (eval = FALSE)
###################################################
## system.time({
## MoranNb <- lapply(vars, function(x) f_bpara(x=eire_ge1[[x]], nsim=nsim, listw=nb_B))
## MoranRb <- lapply(vars, function(x) f_bperm(x=eire_ge1[[x]], nsim=nsim, listw=nb_B))
## Prop_unstdNb  <- lapply(vars, function(x) f_bpara(x=eire_ge1[[x]], nsim=nsim, listw=lw_unstand))
## Prop_unstdRb  <- lapply(vars, function(x) f_bperm(x=eire_ge1[[x]], nsim=nsim, listw=lw_unstand))
## Prop_stdNb  <- lapply(vars, function(x) f_bpara(x=eire_ge1[[x]], nsim=nsim, listw=lw_std))
## Prop_stdRb  <- lapply(vars, function(x) f_bperm(x=eire_ge1[[x]], nsim=nsim, listw=lw_std))
## })


###################################################
### code chunk number 37: CO69.Rnw:626-639 (eval = FALSE)
###################################################
## zzz <- system.time({
## MoranNb <- lapply(vars, function(x) f_bpara(x=eire_ge1[[x]], nsim=nsim, listw=nb_B))
## MoranRb <- lapply(vars, function(x) f_bperm(x=eire_ge1[[x]], nsim=nsim, listw=nb_B))
## Prop_unstdNb  <- lapply(vars, function(x) f_bpara(x=eire_ge1[[x]], nsim=nsim, listw=lw_unstand))
## Prop_unstdRb  <- lapply(vars, function(x) f_bperm(x=eire_ge1[[x]], nsim=nsim, listw=lw_unstand))
## Prop_stdNb  <- lapply(vars, function(x) f_bpara(x=eire_ge1[[x]], nsim=nsim, listw=lw_std))
## Prop_stdRb  <- lapply(vars, function(x) f_bperm(x=eire_ge1[[x]], nsim=nsim, listw=lw_std))
## })
## res <- lapply(c("MoranNb", "MoranRb", "Prop_unstdNb", "Prop_unstdRb", "Prop_stdNb", "Prop_stdRb"), function(x) sapply(get(x), function(y) (y$t0 - mean(y$t))/sd(y$t)))
## res <- t(do.call("rbind", res))
## colnames(res) <- c("MoranNb", "MoranRb", "Prop_unstdNb", "Prop_unstdRb", "Prop_stdNb", "Prop_stdRb")
## rownames(res) <- vars
## save(zzz, res, file="backstore/boot_res.RData")


###################################################
### code chunk number 38: CO69.Rnw:641-644
###################################################
bsfn <- system.file("doc/backstore/boot_res.RData", package="spdep")
load(bsfn)
zzz


###################################################
### code chunk number 39: CO69.Rnw:646-650 (eval = FALSE)
###################################################
## res <- lapply(c("MoranNb", "MoranRb", "Prop_unstdNb", "Prop_unstdRb", "Prop_stdNb", "Prop_stdRb"), function(x) sapply(get(x), function(y) (y$t0 - mean(y$t))/sd(y$t)))
## res <- t(do.call("rbind", res))
## colnames(res) <- c("MoranNb", "MoranRb", "Prop_unstdNb", "Prop_unstdRb", "Prop_stdNb", "Prop_stdRb")
## rownames(res) <- vars


###################################################
### code chunk number 40: CO69.Rnw:664-669
###################################################
print(formatC(res, format="f", digits=4), quote=FALSE)
oores <- ores - res
apply(oores, 2, mad)
alpha_0.05 <- qnorm(0.05, lower.tail=FALSE)
all((res >= alpha_0.05) == (ores >= alpha_0.05))


###################################################
### code chunk number 41: CO69.Rnw:685-693
###################################################
lm_objs <- lapply(vars, function(x) lm(formula(paste(x, "~1")), data=eire_ge1))
system.time({
MoranSad <- lapply(lm_objs, function(x) lm.morantest.sad(x, listw=nb_B))
Prop_unstdSad  <- lapply(lm_objs, function(x) lm.morantest.sad(x, listw=lw_unstand))
Prop_stdSad  <- lapply(lm_objs, function(x) lm.morantest.sad(x, listw=lw_std))
})
res <- sapply(c("MoranSad", "Prop_unstdSad", "Prop_stdSad"), function(x) sapply(get(x), "[[", "statistic"))
rownames(res) <- vars


###################################################
### code chunk number 42: CO69.Rnw:709-713
###################################################
print(formatC(res, format="f", digits=4), quote=FALSE)
oores <- res - ores[,c(1,3,5)]
apply(oores, 2, mad)
all((res >= alpha_0.05) == (ores[,c(1,3,5)] >= alpha_0.05))


###################################################
### code chunk number 43: CO69.Rnw:724-731
###################################################
system.time({ 
MoranEx <- lapply(lm_objs, function(x) lm.morantest.exact(x, listw=nb_B))
Prop_unstdEx  <- lapply(lm_objs, function(x) lm.morantest.exact(x, listw=lw_unstand))
Prop_stdEx  <- lapply(lm_objs, function(x) lm.morantest.exact(x, listw=lw_std))
})
res <- sapply(c("MoranEx", "Prop_unstdEx", "Prop_stdEx"), function(x) sapply(get(x), "[[", "statistic"))
rownames(res) <- vars


###################################################
### code chunk number 44: CO69.Rnw:740-744
###################################################
print(formatC(res, format="f", digits=4), quote=FALSE)
oores <- res - ores[,c(1,3,5)]
apply(oores, 2, mad)
all((res >= alpha_0.05) == (ores[,c(1,3,5)] >= alpha_0.05))


###################################################
### code chunk number 45: CO69.Rnw:756-765
###################################################
vars_scaled <- lapply(vars, function(x) scale(eire_ge1[[x]], scale=FALSE))
nb_W <- nb2listw(lw_unstand$neighbours, style="W")
pre <- spdep:::preAple(0, listw=nb_W)
MoranAPLE <- sapply(vars_scaled, function(x) spdep:::inAple(x, pre))
pre <- spdep:::preAple(0, listw=lw_std, override_similarity_check=TRUE)
Prop_stdAPLE <- sapply(vars_scaled, function(x) spdep:::inAple(x, pre))
res <- cbind(MoranAPLE, Prop_stdAPLE)
colnames(res) <- c("APLE W", "APLE Gstd")
rownames(res) <- vars


###################################################
### code chunk number 46: CO69.Rnw:781-782
###################################################
print(formatC(res, format="f", digits=4), quote=FALSE)


###################################################
### code chunk number 47: CO69.Rnw:788-807
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
opar <- par(mar=c(3,3,1,1)+0.1)
pal <- grey.colors(9, 1, 0.5, 2.2)
oopar <- par(mfrow=c(1,3), mar=c(1,1,3,1)+0.1)
z <- t(listw2mat(nb_B))
brks <- c(0,0.1,1)
image(1:25, 1:25, z[,ncol(z):1], breaks=brks, col=pal[c(1,9)],
 main="Binary", axes=FALSE)
box()
z <- t(listw2mat(lw_unstand))
brks <- c(0,quantile(c(z)[c(z) > 0], seq(0,1,1/8)))
image(1:25, 1:25, z[,ncol(z):1], breaks=brks, col=pal, main="General", axes=FALSE)
box()
z <- t(listw2mat(lw_std))
brks <- c(0,quantile(c(z)[c(z) > 0], seq(0,1,1/8)))
image(1:25, 1:25, z[,ncol(z):1], breaks=brks, col=pal,
 main="Std. general", axes=FALSE)
box()
par(oopar)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 48: CO69.Rnw:816-826
###################################################
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
eire_ge1$nb_B <- sapply(nb_B$weights, sum)
eire_ge1$lw_unstand <- sapply(lw_unstand$weights, sum)
library(lattice)
trellis.par.set(sp.theme())
p1 <- spplot(eire_ge1, c("nb_B"), main="Binary")
p2 <- spplot(eire_ge1, c("lw_unstand"), main="General")
print(p1, split=c(1,1,2,1), more=TRUE)
print(p2, split=c(2,1,2,1), more=FALSE)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 49: CO69.Rnw:843-845
###################################################
options("width"=owidth)
options("warn"=ow)


