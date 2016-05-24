### R code from vignette source 'Anx6.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: Anx6.Rnw:135-139
###################################################
owidth <- getOption("width") # largeur des sorties
options(width=60, continue="+ ","warn"=-1 )
.PngNo <- 0
nom.fich = "./Figures/anx6-bitmap-"


###################################################
### code chunk number 2: bfig (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(rep.ima,nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 7, height = 7, pointsize = 12, bg = "white")


###################################################
### code chunk number 3: bfigps (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""), width = 7, height = 7, pointsize = 12, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 4: bfig1 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(rep.ima,nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 5, height = 2, pointsize = 10, bg = "white")


###################################################
### code chunk number 5: bfigps1 (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""),  width = 5, height =2, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 6: bfig2 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(rep.ima,nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 3.9, height = 3.1, pointsize = 10, bg = "white")


###################################################
### code chunk number 7: bfigps2 (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""), width = 3.9, height = 3.1,   pointsize = 10, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 8: bfig3 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(rep.ima,nom.fich, .PngNo, sep="")
## pdf(file=paste(file,".pdf",sep=""), width = 5.92, height = 6.74, pointsize = 10, bg = "white")


###################################################
### code chunk number 9: bfigps3 (eval = FALSE)
###################################################
## postscript(file=paste(file,".ps",sep=""), width = 5.92, height = 6.74, pointsize = 10, bg = "white",horizontal= FALSE,paper="special")


###################################################
### code chunk number 10: bfig4 (eval = FALSE)
###################################################
## .PngNo <- .PngNo + 1; file = paste(rep.ima,nom.fich, .PngNo, sep="")
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
### code chunk number 14: alpha
###################################################
alpha= seq(.1,.3, by=.1)
arret = seq(10,40,by=10)
n.al= length(alpha) ; n.arret = length(arret)
cumul = matrix(0,nrow=n.al,ncol=n.arret)
rownames(cumul) =as.character(alpha)
colnames(cumul) = as.character(arret)
poids = function(alf,i)
{
# renvoie les poids alpha*(1 - alpha)^j, j=0, i-1
wgh= rep(0,i)
wgh[1]= alf
for(k in 2:i )
{wgh[k] = wgh[k-1]*(1 - alf)}
sum(wgh)
}
for (m in 1:length(alpha))
{
for (n in 1:length(arret))
{
cumul[m,n] = poids(alpha[m],arret[n])
}
}
round(cumul,digits=2)


###################################################
### code chunk number 15: Anx6.Rnw:242-248
###################################################
require(forecast)
require(expsmooth)
require(caschrono)
ets0 = ets(fmsales,model="ANN")
summary(ets0)
str(ets0, width = 60, strict.width = "cut")


###################################################
### code chunk number 16: testbl
###################################################
Box.test.2(residuals(ets0), nlag = c(3,6,9))


###################################################
### code chunk number 17: Anx6.Rnw:289-291
###################################################
(ets0.hw=HoltWinters(fmsales, alpha = NULL, beta = FALSE, 
gamma =FALSE))


