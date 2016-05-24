### R code from vignette source 'Anx7.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: Anx7.Rnw:135-139
###################################################
owidth <- getOption("width") # largeur des sorties
options(width=60, continue="+ ","warn"=-1 )
.PngNo <- 0
nom.fich = "./Figures/annexe_simul-bitmap-"


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
### code chunk number 14: simusar.ini
###################################################
# graine
set.seed(2761)
innov1 = rnorm(290,sd=4.18)
y = arima.sim(list(order = c(12,0,1), ma=-.7, ar=c(rep(0,11),.9)), 
  innov =innov1, n.start =50, n = 240) +50
y.ts = ts(y,frequency=12, start=c(1920,1))
ytr=cbind(y.ts,nottem)
colnames(ytr)=c("serie simulee","temperature")


###################################################
### code chunk number 15: autop
###################################################
require(polynom)
autop=polynomial(c(1,-1/1.4))*polynomial(c(1,-1))*polynomial(c(1,-1/1.9))


###################################################
### code chunk number 16: autop2
###################################################
autop1 = polynomial(c(1,-1/1.4))*polynomial(c(1,-1/1.9))
asim8b = arima.sim(n=60, list(ar = -autop1[-1], 
order=c(2,1,0)))


###################################################
### code chunk number 17: construct3
###################################################
require(dse)
AR = array( autop1,c(length(autop1),1,1))
MA = array(1,c(1,1,1))
mod2 = ARMA(A=AR  ,B = MA)
asim8c= simulate(mod2,sampleT=60, sd=1.5)


