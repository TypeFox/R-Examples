### R code from vignette source 'spp-dist.Rnw'

###################################################
### code chunk number 1: spp-dist.Rnw:1-3
###################################################
options(width=70)
options(continue=" ")


###################################################
### code chunk number 2: spp-dist.Rnw:58-60
###################################################
library(unmarked)
library(raster)


###################################################
### code chunk number 3: spp-dist.Rnw:86-94
###################################################
data(crossbill)
umf <- unmarkedFrameOccu(
    y=as.matrix(crossbill[,c("det991", "det992", "det993")]),
    siteCovs=crossbill[,c("ele", "forest")],
    obsCovs=list(date=crossbill[,c("date991", "date992", "date993")]))
sc <- scale(siteCovs(umf))
siteCovs(umf) <- sc
head(umf)


###################################################
### code chunk number 4: spp-dist.Rnw:110-111
###################################################
(fm.occu <- occu(~date ~ele + I(ele^2) + forest, umf))


###################################################
### code chunk number 5: swiss
###################################################
data(Switzerland)
print(levelplot(elevation ~ x + y, Switzerland, aspect="iso",
                xlab="Easting (m)", ylab="Northing (m)",
                col.regions=terrain.colors(100)))


###################################################
### code chunk number 6: spp-dist.Rnw:137-143
###################################################
library(raster)
elevation <- rasterFromXYZ(Switzerland[,c("x","y","elevation")],
    crs="+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs")

forest <- rasterFromXYZ(Switzerland[,c("x","y","forest")],
    crs="+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs")


###################################################
### code chunk number 7: ef
###################################################
attr(sc, "scaled:center")
attr(sc, "scaled:scale")
ele.s <- (elevation-1189)/640
forest.s <- (forest-34.7)/27.7
ef <- stack(ele.s, forest.s)
names(ef) <- c("ele", "forest")
plot(ef, col=terrain.colors(100))


###################################################
### code chunk number 8: psi
###################################################
(beta <- coef(fm.occu, type="state"))
logit.psi <- beta[1] + beta[2]*ele.s + beta[3]*ele.s^2 + beta[4]*forest.s
psi <- exp(logit.psi) / (1 + exp(logit.psi))
#plot(psi, col=terrain.colors(100))
print(spplot(psi, col.regions=terrain.colors(100)))


###################################################
### code chunk number 9: psi2 (eval = FALSE)
###################################################
## E.psi <- predict(fm.occu, type="state", newdata=ef)
## plot(E.psi, axes=FALSE, col=terrain.colors(100))


###################################################
### code chunk number 10: spp-dist.Rnw:282-294
###################################################
data(issj)
covs <- scale(issj[,c("elevation", "forest", "chaparral")])
area <- pi*300^2 / 10000
jayumf <- unmarkedFrameDS(y=as.matrix(issj[,1:3]),
                          siteCovs=data.frame(covs, area),
                          dist.breaks=c(0,100,200,300),
                          unitsIn="m", survey="point")
head(jayumf)
fm1 <- distsamp(~chaparral ~chaparral + elevation + offset(log(area)),
                jayumf, keyfun="halfnorm", output="abund",
                starts=c(-2.8,1,0,4.5,0))
fm1


###################################################
### code chunk number 11: spp-dist.Rnw:306-312
###################################################
data(issj)
data(cruz)
elev <- rasterFromXYZ(cruz[,c("x","y","elevation")],
     crs="+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
#plot(elev, col=terrain.colors(100))
#points(issj[,c("x","y")], cex=0.5)


###################################################
### code chunk number 12: spp-dist.Rnw:334-350
###################################################
data(cruz)
elev <- rasterFromXYZ(cruz[,c("x","y","elevation")],
     crs="+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
forest <- rasterFromXYZ(cruz[,c("x","y","forest")],
     crs="+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
chap <- rasterFromXYZ(cruz[,c("x","y","chaparral")],
     crs="+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
area.raster <- chap
values(area.raster) <- 300*300/10000 # area of a grid pixel
attr(covs, "scaled:center")
attr(covs, "scaled:scale")
elev.s <- (elev-202)/125
forest.s <- (forest-0.0673)/0.137
chap.s <- (chap-0.270)/0.234
habitat <- stack(elev.s, forest.s, chap.s, area.raster)
names(habitat) <- c("elevation", "forest", "chaparral", "area")


###################################################
### code chunk number 13: issj
###################################################
E <- predict(fm1, type="state", newdata=habitat)
plot(E, axes=FALSE, col=terrain.colors(100))


###################################################
### code chunk number 14: spp-dist.Rnw:374-384
###################################################
cruz2 <- data.frame(cruz[,1:2],
                    chaparral=(cruz$chaparral-0.270)/0.234,
                    elevation=(cruz$elevation-202)/125)
cruz2$E.N <- exp(-2.827 + 0.957*cruz2$chaparral + -0.244*cruz2$elevation)
wireframe(E.N ~ x + y, cruz2,
          shade=TRUE, #shade.colors.palette=terrain.colors(100),
#          drape=TRUE,
          aspect=0.5, colorkey=FALSE,
          screen=list(z=10, x=-10))



###################################################
### code chunk number 15: spp-dist.Rnw:390-391
###################################################
detach(package:raster)


