### R code from vignette source 'Anx10.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: Anx10.Rnw:135-139
###################################################
owidth <- getOption("width") # largeur des sorties
options(width=60, continue="+ ","warn"=-1 )
.PngNo <- 0
nom.fich = "anx10-bitmap-"


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
### code chunk number 12: bfig5 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 7, height = 4, pointsize = 10, bg = "white")


###################################################
### code chunk number 13: bfigps5 (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""),  width = 7, height =4, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 14: zfig2 (eval = FALSE)
###################################################
## dev.null <- dev.off()


###################################################
### code chunk number 15: zfiginclude (eval = FALSE)
###################################################
## cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 16: init.consel
###################################################
require(caschrono)
data(khct)


###################################################
### code chunk number 17: prep.df
###################################################
khct.df<-as.data.frame(window(cbind(khct,time(khct),
(time(khct)-1977)^2), end=c(1983,12)))
colnames(khct.df) <- c("kwh","htdd","cldd","t1","t1.2")


###################################################
### code chunk number 18: prep.df2
###################################################
mod2 = lm(sqrt(kwh) ~ htdd + cldd + t1 + t1.2, data=khct.df)
u = ts(residuals(mod2), start=c(1970,1), frequency=12)


###################################################
### code chunk number 19: prep.conselec2
###################################################
kwh1rc = window(sqrt(khct[,"kwh"]), end=c(1983,12))
xreg1 = khct.df[ ,c("htdd","cldd","t1","t1.2")]
xreg2 =   xreg1[,-4]
(mdarx3c=Arima(kwh1rc,order=c(1,0,0),seasonal=list(order=c(1,0,1)),
  xreg=xreg2))
u.3c=kwh1rc-as.matrix(xreg2)%*%as.matrix(mdarx3c$coef[5:7])-mdarx3c$coef[4]


###################################################
### code chunk number 20: conselec.super (eval = FALSE)
###################################################
## plot.ts(cbind(u,u.3c),plot.type = "single",lty=1:2)
## abline(h=0)


###################################################
### code chunk number 21: conselec0
###################################################
.PngNo <- .PngNo + 1; file = paste(nom.fich, .PngNo, sep="")
pdf(file=paste(file,".pdf",sep=""), width = 7, height = 4, pointsize = 10, bg = "white")
plot.ts(cbind(u,u.3c),plot.type = "single",lty=1:2)
abline(h=0)
dev.null <- dev.off()
postscript(file=paste(file,".ps",sep=""),  width = 7, height =4, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")
dev.null <- dev.off()
cat("\\includegraphics[width=0.9\\textwidth]{", file, "}\n\n", sep="")


###################################################
### code chunk number 22: Anx10.Rnw:322-327
###################################################
vcovar= mdarx3c$var.coef[c('cldd','htdd'),c('cldd','htdd')]
delta0 = mdarx3c$coef['cldd']- 10*mdarx3c$coef['htdd']
sig2 = t(matrix(c(1,-10)))%*%vcovar%*%matrix(c(1,-10))
(t0= delta0/sig2^.5)
(p.val = 2*(1-pnorm(t0)))


