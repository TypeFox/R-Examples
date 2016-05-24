### R code from vignette source 'datasets.Rnw'

###################################################
### code chunk number 1: datasets.Rnw:5-6
###################################################
options(SweaveHooks=list(fig=function() par(mar=c(1,1,1,1))))


###################################################
### code chunk number 2: datasets.Rnw:25-31
###################################################
library(spatstat)
sdate <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Date")
sversion <- read.dcf(file = system.file("DESCRIPTION", package = "spatstat"),
         fields = "Version")
options(useFancyQuotes=FALSE)


###################################################
### code chunk number 3: datasets.Rnw:201-219
###################################################
opa <- par()
## How to set all margins to zero and eliminate all outer spaces
zeromargins <- function() {
  par(
      mar=rep(0,4),
      omd=c(0,1,0,1),
      xaxs="i",
      yaxs="i"
  )
  invisible(NULL)
}
## Set 'mar'
setmargins <- function(...) {
  x <- c(...)
  x <- rep(x, 4)[1:4]
  par(mar=x)
  invisible(NULL)
}


###################################################
### code chunk number 4: datasets.Rnw:228-229 (eval = FALSE)
###################################################
## plot(amacrine)


###################################################
### code chunk number 5: datasets.Rnw:231-233
###################################################
getOption("SweaveHooks")[["fig"]]()
setmargins(0,1,2,0)
plot(amacrine)


###################################################
### code chunk number 6: datasets.Rnw:242-243 (eval = FALSE)
###################################################
## plot(anemones, markscale=1)


###################################################
### code chunk number 7: datasets.Rnw:245-247
###################################################
getOption("SweaveHooks")[["fig"]]()
setmargins(0,0,2,0)
plot(anemones, markscale=1)


###################################################
### code chunk number 8: datasets.Rnw:260-261 (eval = FALSE)
###################################################
## ants.extra$plotit()


###################################################
### code chunk number 9: datasets.Rnw:263-265
###################################################
getOption("SweaveHooks")[["fig"]]()
setmargins(0,0,1,0)
ants.extra$plotit()


###################################################
### code chunk number 10: datasets.Rnw:273-274
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(austates)


###################################################
### code chunk number 11: datasets.Rnw:284-286 (eval = FALSE)
###################################################
## plot(bdspots, equal.scales=TRUE, pch="+", 
##      panel.args=function(i)list(cex=c(0.15, 0.2, 0.7)[i]))


###################################################
### code chunk number 12: datasets.Rnw:288-292
###################################################
getOption("SweaveHooks")[["fig"]]()
zeromargins()
plot(bdspots, equal.scales=TRUE, pch="+", main="",
     mar.panel=0, hsep=1,
     panel.args=function(i)list(cex=c(0.15, 0.2, 0.7)[i]))


###################################################
### code chunk number 13: datasets.Rnw:302-304 (eval = FALSE)
###################################################
## plot(bei.extra$elev, main="Beilschmiedia")
## plot(bei, add=TRUE, pch=16, cex=0.3)


###################################################
### code chunk number 14: datasets.Rnw:306-309
###################################################
getOption("SweaveHooks")[["fig"]]()
setmargins(0,0,2,0)
plot(bei.extra$elev, main="Beilschmiedia")
plot(bei, add=TRUE, pch=16, cex=0.3)


###################################################
### code chunk number 15: datasets.Rnw:312-318
###################################################
getOption("SweaveHooks")[["fig"]]()
M <- persp(bei.extra$elev, 
           theta=-45, phi=18, expand=7,
           border=NA, apron=TRUE, shade=0.3, 
           box=FALSE, visible=TRUE,
           main="")
perspPoints(bei, Z=bei.extra$elev, M=M, pch=16, cex=0.3)


###################################################
### code chunk number 16: datasets.Rnw:327-328
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(betacells)


###################################################
### code chunk number 17: datasets.Rnw:333-334
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(bramblecanes, cols=1:3)


###################################################
### code chunk number 18: datasets.Rnw:337-338
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(split(bramblecanes))


###################################################
### code chunk number 19: datasets.Rnw:348-349
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(bronzefilter,markscale=2)


###################################################
### code chunk number 20: datasets.Rnw:358-359
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(cells)


###################################################
### code chunk number 21: datasets.Rnw:368-371
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(chicago, main="Chicago Street Crimes", col="grey",
     cols=c("red", "blue", "black", "blue", "red", "blue", "blue"),
     chars=c(16,2,22,17,24,15,6), leg.side="left", show.window=FALSE)


###################################################
### code chunk number 22: datasets.Rnw:381-382
###################################################
getOption("SweaveHooks")[["fig"]]()
chorley.extra$plotit()


###################################################
### code chunk number 23: datasets.Rnw:398-400
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(clmfires, which.marks="cause", cols=2:5, cex=0.25,
     main="Castilla-La Mancha forest fires")


###################################################
### code chunk number 24: datasets.Rnw:410-411
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(clmfires.extra$clmcov200, main="Covariates for forest fires")


###################################################
### code chunk number 25: datasets.Rnw:422-424
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(copper$Points, main="Copper")
plot(copper$Lines, add=TRUE)


###################################################
### code chunk number 26: datasets.Rnw:431-433
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(demohyper, quote({ plot(Image, main=""); plot(Points, add=TRUE) }),
      parargs=list(mar=rep(1,4)))


###################################################
### code chunk number 27: datasets.Rnw:440-441
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(demopat)


###################################################
### code chunk number 28: datasets.Rnw:455-456
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(dendrite, leg.side="bottom", main="", cex=0.75, cols=2:4)


###################################################
### code chunk number 29: datasets.Rnw:464-465
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(finpines, main="Finnish pines")


###################################################
### code chunk number 30: datasets.Rnw:478-482
###################################################
getOption("SweaveHooks")[["fig"]]()
wildM1 <- with(flu, virustype == "wt" & stain == "M2-M1")
plot(flu[wildM1, 1, drop=TRUE],
     main=c("flu data", "wild type virus, M2-M1 stain"),
     chars=c(16,3), cex=0.4, cols=2:3)


###################################################
### code chunk number 31: datasets.Rnw:490-491
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(gordon, main="People in Gordon Square", pch=16)


###################################################
### code chunk number 32: datasets.Rnw:506-507
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(gorillas, which.marks=1, chars=c(1,3), cols=2:3, main="Gorilla nest sites")


###################################################
### code chunk number 33: datasets.Rnw:511-512 (eval = FALSE)
###################################################
## system.file("rawdata/gorillas/vegetation.asc", package="spatstat")


###################################################
### code chunk number 34: datasets.Rnw:521-522
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(hamster, cols=c(2,4))


###################################################
### code chunk number 35: datasets.Rnw:532-533
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(heather)


###################################################
### code chunk number 36: datasets.Rnw:543-544
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(humberside)


###################################################
### code chunk number 37: datasets.Rnw:556-557
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(hyytiala, cols=2:5)


###################################################
### code chunk number 38: datasets.Rnw:566-567
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(japanesepines)


###################################################
### code chunk number 39: datasets.Rnw:576-577
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(lansing)


###################################################
### code chunk number 40: datasets.Rnw:580-581
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(split(lansing))


###################################################
### code chunk number 41: datasets.Rnw:588-589
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(longleaf)


###################################################
### code chunk number 42: datasets.Rnw:598-600
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(mucosa, chars=c(1,3), cols=c("red", "green"))
plot(mucosa.subwin, add=TRUE, lty=3)


###################################################
### code chunk number 43: datasets.Rnw:614-617
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(murchison$greenstone, main="Murchison data", col="lightgreen")
plot(murchison$gold, add=TRUE, pch=3, col="blue")
plot(murchison$faults, add=TRUE, col="red")


###################################################
### code chunk number 44: datasets.Rnw:625-626
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(nbfires, use.marks=FALSE, pch=".")


###################################################
### code chunk number 45: datasets.Rnw:629-630
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(split(nbfires), use.marks=FALSE, chars=".")


###################################################
### code chunk number 46: datasets.Rnw:633-638
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar=c(0,0,2,0))
plot(split(nbfires)$"2000", which.marks="fire.type",
     main=c("New Brunswick fires 2000", "by fire type"),
     cols=c("blue", "green", "red", "cyan"),
     leg.side="left")


###################################################
### code chunk number 47: datasets.Rnw:646-648
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(nztrees)
plot(trim.rectangle(as.owin(nztrees), c(0,5), 0), add=TRUE, lty=3)


###################################################
### code chunk number 48: datasets.Rnw:661-662
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(osteo[1:10,], main.panel="", pch=21, bg='white')


###################################################
### code chunk number 49: datasets.Rnw:668-669 (eval = FALSE)
###################################################
## system.file("rawdata/osteo/osteo36.txt", package="spatstat")


###################################################
### code chunk number 50: datasets.Rnw:678-679
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(paracou, cols=2:3, chars=c(16,3))


###################################################
### code chunk number 51: datasets.Rnw:687-688
###################################################
getOption("SweaveHooks")[["fig"]]()
ponderosa.extra$plotit()


###################################################
### code chunk number 52: datasets.Rnw:699-702
###################################################
getOption("SweaveHooks")[["fig"]]()
pyr <- pyramidal
pyr$grp <- abbreviate(pyramidal$group, minlength=7)
plot(pyr, quote(plot(Neurons, pch=16, main=grp)), main="Pyramidal Neurons")


###################################################
### code chunk number 53: datasets.Rnw:722-724
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(redwood)
plot(redwood3, add=TRUE, pch=20)


###################################################
### code chunk number 54: datasets.Rnw:727-728
###################################################
getOption("SweaveHooks")[["fig"]]()
redwoodfull.extra$plotit()


###################################################
### code chunk number 55: datasets.Rnw:742-744
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(as.listof(residualspaper[c("Fig1", "Fig4a", "Fig4b", "Fig4c")]), 
     main="")


###################################################
### code chunk number 56: datasets.Rnw:752-753
###################################################
getOption("SweaveHooks")[["fig"]]()
shapley.extra$plotit(main="Shapley")


###################################################
### code chunk number 57: datasets.Rnw:760-761
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(simdat)


###################################################
### code chunk number 58: datasets.Rnw:769-770
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(spiders, pch=16, show.window=FALSE)


###################################################
### code chunk number 59: datasets.Rnw:777-780
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(sporophores, chars=c(16,1,2), cex=0.6)
points(0,0,pch=16, cex=2)
text(15,8,"Tree", cex=0.75)


###################################################
### code chunk number 60: datasets.Rnw:789-790
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(spruces, maxsize=min(nndist(spruces)))


###################################################
### code chunk number 61: datasets.Rnw:799-800
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(swedishpines)


###################################################
### code chunk number 62: datasets.Rnw:809-810
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(urkiola, cex=0.5, cols=2:3)


###################################################
### code chunk number 63: datasets.Rnw:817-819
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar=c(0,0,2,0))
plot(waka, markscale=0.04, main=c("Waka national park", "tree diameters"))


###################################################
### code chunk number 64: datasets.Rnw:826-830
###################################################
getOption("SweaveHooks")[["fig"]]()
v <- rotate(vesicles, pi/2)
ve <- lapply(vesicles.extra, rotate, pi/2)
plot(v, main="Vesicles")
plot(ve$activezone, add=TRUE, lwd=3)


###################################################
### code chunk number 65: datasets.Rnw:855-856 (eval = FALSE)
###################################################
## system.file("rawdata/vesicles/mitochondria.txt", package="spatstat")


###################################################
### code chunk number 66: datasets.Rnw:864-865
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(waterstriders)


