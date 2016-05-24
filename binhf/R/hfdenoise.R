`hfdenoise` <-
function(n=256, proportion=P2, binsize=1, thrule="ebayesthresh", van=8, fam="DaubLeAsymm", 
pl=3,prior="laplace",vscale="independent", plotstep=FALSE,truncate=FALSE,...){

   
    x <- (1:n)/n
    y <- proportion(x,...)

J<-log2(n)
    
    binary<-NULL
 
    for (j in 1:n) {
           # binary[j] <- rbinom(1, binsize, y[j]/binsize)
            binary[j] <- rbinom(1, binsize, y[j])
        }

if (plotstep==TRUE){
	plot(x, y, main="Proportion")
	scan()
	plot(x, binary, main="Binary")
	scan()
	}


        hfx1 <- binhf.wd(binary, binsize)
        hfx <- hfx1$transformed
        hfxa <- ansc(binary, binsize)
        hfxf <- free(binary,binsize)

if (plotstep==TRUE)	{
plot(x, hfx, main="HF transformed binary")
scan()
plot(x, hfxa, main="Ans transformed binary")
scan()
}


tmpwd <- wd(hfx, filter.number = van, family = fam)
tmpwda<- wd(hfxa, filter.number = van, family = fam)
tmpwdf<- wd(hfxf, filter.number = van, family = fam)

if (plotstep==TRUE)	{
plot(tmpwd, main="HF transformed wd")
scan()
plot(tmpwda, main="Ans transformed wd")
scan()
}


if (thrule == "sureshrink") {
tmpwd.thresh <- threshold.wd(tmpwd, levels = pl:(tmpwd$nlevels - 
            1), policy = "sure", type = "soft", dev = madmad)
tmpwd.thresha <- threshold.wd(tmpwda, levels = pl:(tmpwda$nlevels - 
            1), policy = "sure", type = "soft", dev = madmad)
tmpwd.threshf <- threshold.wd(tmpwdf, levels = pl:(tmpwdf$nlevels - 
            1), policy = "sure", type = "soft", dev = madmad)

}

else {
tmpwd.thresh <- ebayesthresh.wavelet(tmpwd, prior = prior, 
            a = NA, smooth.levels = (tmpwd$nlevels - pl), vscale=vscale)
tmpwd.thresha <- ebayesthresh.wavelet(tmpwda, prior = prior, 
            a = NA, smooth.levels = (tmpwda$nlevels - pl), vscale=vscale)
tmpwd.threshf <- ebayesthresh.wavelet(tmpwdf, prior = prior, 
            a = NA, smooth.levels = (tmpwdf$nlevels - pl), vscale=vscale)
}

#
# Inversion
#
if (thrule != "complex")	{    
if (plotstep==TRUE)	{
	plot(tmpwd.thresh, main="HF transformed thresholded wd")
	scan()
	plot(tmpwd.thresha, main="Ans transformed thresholded wd")
	scan()
	}
	tmp.wr <- wr(tmpwd.thresh)
	tmp.wra <- wr(tmpwd.thresha)
	tmp.wrf <- wr(tmpwd.threshf)
	}


fhat <- tmp.wr
fhat <- list(transformed = fhat, cnew = hfx1$cnew)
fhat <- invbinhf.wd(fhat, binsize)/binsize

fhata<- tmp.wra
fhata<- invansc(fhata, binsize)/binsize


fhatf<- tmp.wrf
fhatf <- freeinv(fhatf,binsize)/binsize

if(truncate==TRUE){
fhat[fhat<0]<-0
fhata[fhata<0]<-0
fhatf[fhatf<0]<-0

fhat[fhat>1]<-1
fhata[fhata>1]<-1
}


if (plotstep==TRUE)	{

plot(x, fhat, type="l", main="HF estimate")
scan()
plot(x, fhata, type="l", main="Ansc estimate")
}


l <- list(x=x,truep=y, fhat=fhat, fhata=fhata, fhatf=fhatf,bbwd=tmpwd, awd=tmpwda,b=binary,bb=hfx1,thr=tmpwd.thresh,tmp=tmp.wr)
l

}

