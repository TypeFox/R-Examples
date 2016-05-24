### R code from vignette source 'distsamp.Rnw'

###################################################
### code chunk number 1: distsamp.Rnw:1-3
###################################################
options(width=70)
options(continue=" ")


###################################################
### code chunk number 2: distsamp.Rnw:105-113
###################################################
library(unmarked)
dists <- read.csv(system.file("csv", "distdata.csv", package="unmarked"))
head(dists, 10)
levels(dists$transect) <- c(levels(dists$transect), "g")
levels(dists$transect)
yDat <- formatDistData(dists, distCol="distance",
	transectNameCol="transect", dist.breaks=c(0, 5, 10, 15, 20))
yDat


###################################################
### code chunk number 3: distsamp.Rnw:124-126
###################################################
(covs <- data.frame(canopyHt = c(5, 8, 3, 2, 4, 7, 5),
	habitat = c('A','A','A','A','B','B','B'), row.names=letters[1:7]))


###################################################
### code chunk number 4: distsamp.Rnw:138-141
###################################################
umf <- unmarkedFrameDS(y=as.matrix(yDat), siteCovs=covs, survey="line",
	dist.breaks=c(0, 5, 10, 15, 20), tlength=rep(100, 7),
	unitsIn="m")


###################################################
### code chunk number 5: umfhist
###################################################
summary(umf)
hist(umf, xlab="distance (m)", main="", cex.lab=0.8, cex.axis=0.8)


###################################################
### code chunk number 6: distsamp.Rnw:187-192
###################################################
hn_Null <- distsamp(~1~1, umf)
hn_Null <- distsamp(~1~1, umf, keyfun="halfnorm", output="density",
	unitsOut="ha")
haz_Null <- distsamp(~1~1, umf, keyfun="hazard")
hn_Hab.Ht <- distsamp(~canopyHt ~habitat, umf)


###################################################
### code chunk number 7: distsamp.Rnw:207-208
###################################################
haz_Null


###################################################
### code chunk number 8: distsamp.Rnw:216-221
###################################################
names(haz_Null)
backTransform(haz_Null, type="state")
backTransform(haz_Null, type="det")
backTransform(haz_Null, type="scale")
backTransform(linearComb(hn_Hab.Ht['det'], c(1, 5)))


###################################################
### code chunk number 9: distsamp.Rnw:240-244
###################################################
site.level.density <- predict(hn_Hab.Ht, type="state")$Predicted
plotArea.inHectares <- 100 * 40 / 10000
site.level.abundance <- site.level.density * plotArea.inHectares
(N.hat <- sum(site.level.abundance))


###################################################
### code chunk number 10: distsamp.Rnw:252-260
###################################################
getN.hat <- function(fit) {
    d <- predict(fit, type="state")$Predicted
    a <- d * (100 * 40 / 10000)
    N.hat <- c(N.hat = sum(a))
    return(N.hat)
    }
pb <- parboot(hn_Hab.Ht, statistic=getN.hat, nsim=25)
pb


###################################################
### code chunk number 11: distsamp.Rnw:283-287
###################################################
head(habConstant <- data.frame(canopyHt = seq(2, 8, length=20),
	habitat=factor("A", levels=c("A", "B"))))
(htConstant <- data.frame(canopyHt = 5,
	habitat=factor(c("A", "B"))))


###################################################
### code chunk number 12: distsamp.Rnw:293-297
###################################################
(Elambda <- predict(hn_Hab.Ht, type="state", newdata=htConstant,
	appendData=TRUE))
head(Esigma <- predict(hn_Hab.Ht, type="det", newdata=habConstant,
	appendData=TRUE))


###################################################
### code chunk number 13: predplot
###################################################
par(mfrow=c(1, 2))
with(Elambda, {
	x <- barplot(Predicted, names=habitat, xlab="Habitat",
		ylab="Density (animals / ha)", ylim=c(0, 25), cex.names=0.7,
        cex.lab=0.7, cex.axis=0.7)
	arrows(x, Predicted, x, Predicted+SE, code=3, angle=90, length=0.05)
	box()
	})
with(Esigma, {
	plot(canopyHt, Predicted, type="l", xlab="Canopy height",
		ylab=expression(paste("Half-normal scale parameter (", sigma, ")",
        sep="")), ylim=c(0, 20), cex=0.7, cex.lab=0.7,
    cex.axis=0.7)
	lines(canopyHt, Predicted+SE, lty=2)
	lines(canopyHt, Predicted-SE, lty=2)
	})


###################################################
### code chunk number 14: detplot
###################################################
par(mfrow=c(1, 2))
plot(function(x) gxhn(x, sigma=10.8), 0, 20, xlab="Distance (m)",
	ylab="Detection prob. at 2m canopy ht.", cex.lab=0.7,
    cex.axis=0.7, las=1)
hist(hn_Null, xlab="Distance (m)", ylab="Probability density", main="",
    ylim=c(0, 0.1), cex.lab=0.7, cex.axis=0.7, las=1)
plot(function(x) dxhn(x, sigma=10.8), 0, 20, add=TRUE, col="blue")
plot(function(x) dxhn(x, sigma=9.9), 0, 20, add=TRUE, col="green")
legend('topright', c("Canopy ht. = 2m", "Null", "Canopy ht. = 8m"),
	col=c("blue", "black", "green"), lty=1, cex=0.4)


