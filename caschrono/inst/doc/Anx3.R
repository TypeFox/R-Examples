### R code from vignette source 'Anx3.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: Anx3.Rnw:137-142
###################################################
owidth <- getOption("width") # largeur des sorties
options(width=60, continue="+ ","warn"=-1 )
.PngNo <- 0
#rep.ima = getwd()
nom.fich = "anx3-bitmap-"


###################################################
### code chunk number 2: bfig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 7, height = 7, pointsize = 12, bg = "white")


###################################################
### code chunk number 3: bfigps (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""), width = 7, height = 7, pointsize = 12, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 4: bfig1 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 5, height = 2, pointsize = 10, bg = "white")


###################################################
### code chunk number 5: bfigps1 (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""),  width = 5, height =2, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 6: bfig2 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 3.9, height = 3.1, pointsize = 10, bg = "white")


###################################################
### code chunk number 7: bfigps2 (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""), width = 3.9, height = 3.1,   pointsize = 10, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 8: bfig3 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 5.92, height = 6.74, pointsize = 10, bg = "white")


###################################################
### code chunk number 9: bfigps3 (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""), width = 5.92, height = 6.74, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 10: bfig4 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 6, height = 6, pointsize = 10, bg = "white")


###################################################
### code chunk number 11: bfigps4 (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""), width = 6, height = 6, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 12: zfig2 (eval = FALSE)
###################################################
## dev.null <- dev.off()


###################################################
### code chunk number 13: zfiginclude (eval = FALSE)
###################################################
## cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 14: sim.nor
###################################################
require(e1071)
set.seed(5923)
x= rnorm(1000)
xr = exp(x)
c(skewness(x),kurtosis(x))
c(skewness(xr),kurtosis(xr))


###################################################
### code chunk number 15: plot.nor (eval = FALSE)
###################################################
## op=par(mfrow=c(1,2),mar=c(4,3,2,0.5),lwd=2 )
## plot(density(x),main="",ylab="density",xlab="x normale")
## plot(density(xr),main="",ylab="density",xlab="xr log-normale")
## par(op)


###################################################
### code chunk number 16: Anx3.Rnw:459-466
###################################################
.PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
pdf(file=paste(file,".pdf",sep=""), width = 6, height = 6, pointsize = 10, bg = "white")
op=par(mfrow=c(1,2),mar=c(4,3,2,0.5),lwd=2 )
plot(density(x),main="",ylab="density",xlab="x normale")
plot(density(xr),main="",ylab="density",xlab="xr log-normale")
par(op)
dev.null <- dev.off()
postscript(file=paste(file,".ps",sep=""), width = 6, height = 6, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")
op=par(mfrow=c(1,2),mar=c(4,3,2,0.5),lwd=2 )
plot(density(x),main="",ylab="density",xlab="x normale")
plot(density(xr),main="",ylab="density",xlab="xr log-normale")
par(op)
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


