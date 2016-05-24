### R code from vignette source 'SpatialFiltering.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: SpatialFiltering.Rnw:42-47
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
### code chunk number 8: SpatialFiltering.Rnw:112-120
###################################################
library(maptools)
library(spdep)
owd <- getwd()
setwd(system.file("etc/shapes", package="spdep"))
NY8 <- readShapeSpatial("NY8_utm18")
setwd(system.file("etc/weights", package="spdep"))
NY_nb <- read.gal("NY_nb.gal", region.id=row.names(NY8))
setwd(owd)


###################################################
### code chunk number 9: SpatialFiltering.Rnw:126-131
###################################################
nySFE <- SpatialFiltering(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, nb=NY_nb, style="W", verbose=FALSE)
nylmSFE <- lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME+fitted(nySFE), data=NY8)
summary(nylmSFE)
nylm <- lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8)
anova(nylm, nylmSFE)


###################################################
### code chunk number 10: SpatialFiltering.Rnw:155-157
###################################################
NYlistwW <- nb2listw(NY_nb, style = "W")
set.seed(111)


###################################################
### code chunk number 11: SpatialFiltering.Rnw:158-159 (eval = FALSE)
###################################################
## nyME <- ME(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, offset=log(POP8), family="poisson", listw=NYlistwW, alpha=0.5)


###################################################
### code chunk number 12: SpatialFiltering.Rnw:160-162
###################################################
bsfn <- system.file("doc/backstore/nyME_res.RData", package="spdep")
load(bsfn)


###################################################
### code chunk number 13: SpatialFiltering.Rnw:164-167
###################################################
nyME
NY8$eigen_24 <- fitted(nyME)[,1]
NY8$eigen_223 <- fitted(nyME)[,2]


###################################################
### code chunk number 14: SpatialFiltering.Rnw:174-182
###################################################
.iwidth <- 6
.iheight <- 4
.ipointsize <- 10
.PngNo <- .PngNo + 1; file <- paste("Fig-bitmap-", .PngNo, ".pdf", sep="")
pdf(file=file, width = 6.5, height = 3.5, pointsize = 12, bg = "white")
opar <- par(mar=c(3,3,1,1)+0.1)
library(RColorBrewer)
#gry <- brewer.pal(9, "Greys")[-1]
spplot(NY8, c("eigen_24", "eigen_223"), col.regions=grey.colors(6, 0.95, 0.55, 2.2), cuts=5)
par(opar)
dev.null <- dev.off()
cat("\\includegraphics[width=0.95\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 15: SpatialFiltering.Rnw:190-194
###################################################
nyglmME <- glm(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8))+fitted(nyME), data=NY8, family="poisson")
summary(nyglmME)
nyGLMp <- glm(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8)), data=NY8,family="poisson")
anova(nyGLMp, nyglmME, test="Chisq")


