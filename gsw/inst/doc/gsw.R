### R code from vignette source 'gsw.Rnw'

###################################################
### code chunk number 1: gsw.Rnw:137-138
###################################################
options(keep.source=TRUE, width=60, prompt=' ', continue=' ', oceEOS="unesco")


###################################################
### code chunk number 2: gsw.Rnw:147-149
###################################################
library(gsw)
SA <- gsw_SA_from_SP(SP=35, p=100, longitude=188, latitude=4)


###################################################
### code chunk number 3: gsw.Rnw:153-154
###################################################
CT <- gsw_CT_from_t(SA=SA, t=10, p=100)


###################################################
### code chunk number 4: gsw.Rnw:182-183 (eval = FALSE)
###################################################
## pdf('TS_unesco.pdf', pointsize=18)


###################################################
### code chunk number 5: gsw.Rnw:185-191 (eval = FALSE)
###################################################
## library(oce)
## data(section)
## ctd <- section[["station", 100]]
## Slim <- c(34.8, 37.0)
## Tlim <- c(0, 25)
## plotTS(ctd, Slim=Slim, Tlim=Tlim, eos="unesco")


###################################################
### code chunk number 6: gsw.Rnw:193-195 (eval = FALSE)
###################################################
## dev.off()
## pdf('TS_gsw.pdf', pointsize=18)


###################################################
### code chunk number 7: gsw.Rnw:199-200 (eval = FALSE)
###################################################
## plotTS(ctd, Slim=Slim, Tlim=Tlim, eos="gsw")


###################################################
### code chunk number 8: gsw.Rnw:202-203 (eval = FALSE)
###################################################
## dev.off()


###################################################
### code chunk number 9: gsw.Rnw:224-226 (eval = FALSE)
###################################################
## pdf('temperature_comparison.pdf', height=5, pointsize=18)
## par(mar=c(3.2, 3, 1, 1/2), mgp=c(2, 0.85, 0))


###################################################
### code chunk number 10: gsw.Rnw:228-229 (eval = FALSE)
###################################################
## hist(section[["theta"]] / section[["CT"]], main="")


###################################################
### code chunk number 11: gsw.Rnw:233-236 (eval = FALSE)
###################################################
## dev.off()
## pdf('salinity_comparison.pdf', height=5, pointsize=18)
## par(mar=c(3.2, 3, 1, 1/2), mgp=c(2, 0.85, 0))


###################################################
### code chunk number 12: gsw.Rnw:238-239 (eval = FALSE)
###################################################
## hist(section[["salinity"]] / section[["SA"]], main="")


###################################################
### code chunk number 13: gsw.Rnw:241-242 (eval = FALSE)
###################################################
## dev.off()


###################################################
### code chunk number 14: gsw.Rnw:262-263 (eval = FALSE)
###################################################
## png('SSS_%d.png', width=7, height=4, unit="in", res=150, pointsize=14)


###################################################
### code chunk number 15: gsw.Rnw:265-276 (eval = FALSE)
###################################################
## data("levitus", package="ocedata")
## SSS <- levitus$SSS
## dim <- dim(SSS)
## ll <- expand.grid(lon=levitus$longitude, lat=levitus$latitude)
## SA <- gsw_SA_from_SP(levitus$SSS, 0, ll$lon, ll$lat)
## imagep(levitus$longitude, levitus$latitude, levitus$SSS, col=oceColorsJet)
## title("Surface SP")
## per <- 100 * (1 - levitus$SSS / SA)
## imagep(levitus$longitude, levitus$latitude, per, col=oceColorsJet,
##        zlim=quantile(per, c(0.001, 0.999), na.rm=TRUE))
## title("Surface SA-SP, percent")


###################################################
### code chunk number 16: gsw.Rnw:278-279 (eval = FALSE)
###################################################
## dev.off()


