### R code from vignette source 'presentation_geoxp.rnw'

###################################################
### code chunk number 1: presentation_geoxp.rnw:31-36
###################################################
owidth <- getOption("width")
options("width"=70)
ow <- getOption("warn")
options("warn"=-1)
.PngNo <- 0


###################################################
### code chunk number 2: bfigl (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1;
## file <- paste("Fig", .PngNo, ".pdf", sep="")
## file2 <- paste("Fig", .PngNo, sep="")
## pdf(file=file, width = 8, height = 8, onefile = FALSE, paper = "special")


###################################################
### code chunk number 3: zfigla (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.55\\textwidth]{", file2, "}\n\n", sep="")
## #cat("\\includegraphics[bb = 72 500 288 750]{", file, "}\n\n", sep="")


###################################################
### code chunk number 4: zfiglb (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.55\\textwidth]{", file2, "}\n\n", sep="")
## #cat("\\includegraphics[bb = 72 500 288 750]{", file, "}\n\n", sep="")


###################################################
### code chunk number 5: zfiglc (eval = FALSE)
###################################################
## dev.null <- dev.off()
## cat("\\includegraphics[width=0.9\\textwidth]{", file2, "}\n\n", sep="")
## #cat("\\includegraphics[bb = 72 500 288 750]{", file, "}\n\n", sep="")


###################################################
### code chunk number 6: presentation_geoxp.rnw:67-68
###################################################
library(GeoXp)


###################################################
### code chunk number 7: presentation_geoxp.rnw:110-113
###################################################
data(immob)
class(immob)
names(immob)


###################################################
### code chunk number 8: presentation_geoxp.rnw:125-127
###################################################
immob.sp = SpatialPoints(cbind(immob$longitude,immob$latitude))
class(immob.sp)


###################################################
### code chunk number 9: presentation_geoxp.rnw:135-137
###################################################
immob.spdf = SpatialPointsDataFrame(immob.sp, immob)
class(immob.spdf )


###################################################
### code chunk number 10: presentation_geoxp.rnw:149-150 (eval = FALSE)
###################################################
## histomap(immob.spdf,"prix.vente")


###################################################
### code chunk number 11: presentation_geoxp.rnw:228-231
###################################################
require("maptools")
midiP <- readShapePoly(system.file("shapes/region.shp", package="GeoXp")[1])
cont_midiP<-spdf2list(midiP[-c(22,23),])$poly


###################################################
### code chunk number 12: presentation_geoxp.rnw:238-239
###################################################
criteria <- (immob$latitude>mean(immob$latitude))


###################################################
### code chunk number 13: presentation_geoxp.rnw:250-255 (eval = FALSE)
###################################################
## histomap(immob.spdf, 7, nbcol=15, type = "percent",
## names.attr=names(immob), criteria=criteria, carte=cont_midiP,
## identify=TRUE, cex.lab=0.5, pch=12, col="pink",
## xlab="variation price", ylab="percent", axes=TRUE, lablong="x",
## lablat="y")


###################################################
### code chunk number 14: presentation_geoxp.rnw:280-281
###################################################
last.select<-c(12,18,24,31,32,37,39,42,49,67,73,74,79,81,84)


###################################################
### code chunk number 15: presentation_geoxp.rnw:286-287
###################################################
print("Results have been saved in last.select object")


###################################################
### code chunk number 16: presentation_geoxp.rnw:289-290
###################################################
last.select


###################################################
### code chunk number 17: presentation_geoxp.rnw:363-366 (eval = FALSE)
###################################################
## dbledensitymap(immob.spdf,c("prix.vente","prix.location"),
## xlab=c("selling price","rending price"),identify=TRUE, cex.lab=0.5,
## carte=cont_midiP)


###################################################
### code chunk number 18: presentation_geoxp.rnw:422-424
###################################################
W.nb<-tri2nb(cbind(immob$longitude,immob$latitude))
class(W.nb)


###################################################
### code chunk number 19: presentation_geoxp.rnw:434-437
###################################################
W2.matrix<-makeneighborsw(cbind(immob$longitude,immob$latitude),method="both",m=5,d=175000)
W2.nb<-mat2listw(W2.matrix)$neighbours
class(W2.nb)


###################################################
### code chunk number 20: presentation_geoxp.rnw:450-452 (eval = FALSE)
###################################################
## neighbourmap(immob.spdf,"prix.vente", W.nb, identify=TRUE, cex.lab=0.5,
## carte=cont_midiP)


###################################################
### code chunk number 21: presentation_geoxp.rnw:474-475
###################################################
last.select<-matrix(c(3,17,3,38,3,39,3,40,3,63,3,70,3,90,85,53,85,58,85,65,85,66,85,73),12,2,byrow=TRUE)


###################################################
### code chunk number 22: presentation_geoxp.rnw:480-481
###################################################
print("Results have been saved in last.select object")


###################################################
### code chunk number 23: presentation_geoxp.rnw:483-484
###################################################
last.select


