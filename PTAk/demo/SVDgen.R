#########
 demo.SVDgen <- function(ima=NULL,snr=3,nTens=7,openX11s=FALSE)
#########
  {
# try denoising  an image with svd and with svdsmoo
# ima is an array image more than 25 25
#

library(modreg)

#  how was created timage12
# suj1.12 <- f.get32(name=paste("/usr/people/didier/DATA/Jane/suj",1:12,".img",
#    sep=""),numslice=40,resx=91,resy=109)
# m12 <- apply(suj1.12,c(1,2),mean)
# sd12 <- 1/sqrt(12)*apply(suj1.12,c(1,2),sd)
# mask40 <-f.getmask()
# mask40 <- mask40[,,40]
# timage12 <- Thresh(m12/sd12,mask40,.9,0)
# d.ori <- timage12[,,1]
#

if(is.null(ima)){
        data(timage12)
 # gives d.ori
 # a 12 subjects t-stat  image for a group fmri analysis verbal study

       resx <- 91
       resy <- 109
       }
 else { timage12 <- ima
        resx <- dim(timage12)[1]
        resy <- dim(timage12)[2]
        }
 ## from there you can put any other image in timage12 (with new resx resy)
 # e.g. timage12 <- f.read.analyse.slice("yours.img",1,1)  # uses J Marchini R-package AnalyzeReadR #

if(openX11s)X11(width=4,height=4)
filled.contour(1:resx,1:resy,timage12,main="ori")
 # adding noise

bsr <- range(as.vector(timage12))
bs <- bsr[2]-bsr[1]

cat("----","\n")
cat(" range of the image = ",bs,"\n", "noise added with sd =  range /snr = ",bs/snr,"\n")

    bruitim <- array(rnorm(resx*resy,mean=0,sd=bs/snr),c(resx,resy))
    # too much noise and you will have difficulties to converge with this algorithm

  timage12 <- timage12+bruitim

d.svd <- SVDgen(timage12,nomb=min(25,resx,resy))

  debtime <- proc.time()

  smoofun <- rep(list(Susan1D),7)
  smoofun[[8]] <- NA
d.svdo <- SVDgen(timage12,nomb=min(25,resx,resy),smoothing=TRUE,smoo =list(smoofun))

 duree <- (proc.time()-debtime)[3]

cat("duree ", duree ,"s","\n")

  # took   about  2 mn on my machine
summary(d.svd,testvar=0.5)
summary(d.svdo,testvar=0.5)
if(openX11s)X11(width=6,height=4)
par(mfrow=c(1,2))
plot(d.svd,scree=TRUE,RiskJack=0,nbvs=1:12,main="svd")
plot(d.svdo,scree=TRUE,nbvs=1:12,RiskJack=0,main="svdo")

d.svd.re <- REBUILD(d.svd,nTens=1:nTens,testvar=0)
d.svdo.re <- REBUILD(d.svdo,nTens=1:nTens,testvar=0)

if(openX11s)X11(width=9,height=4)
par(mfrow=c(1,3))
zlimi <- bsr # range(c(range(timage12),range(d.svd.re),range(d.svdo.re)))

image(1:resx,1:resy,timage12,col=cm.colors(20),zlim=zlimi)
image(1:resx,1:resy,d.svd.re,col=cm.colors(20),zlim=zlimi)
image(1:resx,1:resy,d.svdo.re,col=cm.colors(20),zlim=zlimi)
par(mfrow=c(1,1))
if(openX11s)X11(width=4,height=4)
filled.contour(1:resx,1:resy,timage12,main=paste("ori+bruit (snr=",snr,")"))
if(openX11s)X11(width=4,height=4)
filled.contour(1:resx,1:resy,d.svd.re,main=paste("svd",nTens))
if(openX11s)X11(width=4,height=4)
filled.contour(1:resx,1:resy,d.svdo.re,main=paste("svdsmooth",nTens))
}
#########
demo.SVDgen()
cat("\n","\n","args(demo.SVDgen)","\n")
print(args(demo.SVDgen))
