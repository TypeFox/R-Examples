### R code from vignette source 'StatDA.Rnw'

###################################################
### code chunk number 1: fig1plot
###################################################
library(StatDA)
data(chorizon)
Ba=chorizon[,"Ba"]
n=length(Ba)
par(mfcol=c(2,2),mar=c(4,4,2,2))
edaplotlog(Ba,H.freq=F,box=T,H.breaks=30,S.pch=3,S.cex=0.5,D.lwd=1.5,P.log=F,
  P.main="",P.xlab="Ba [mg/kg]",P.ylab="Density",B.pch=3,B.cex=0.5,B.log=TRUE)
edaplot(log10(Ba),H.freq=F,box=T,S.pch=3,S.cex=0.5,D.lwd=1.5,P.ylab="Density",
   P.log=T,P.logfine=c(5,10),P.main="",P.xlab="Ba [mg/kg]",B.pch=3,B.cex=0.5)
plot(sort(Ba),((1:n)-0.5)/n,pch=3,cex=0.8,
main="",xlab="Ba [mg/kg]",ylab="Probability",cex.lab=1,cex.lab=1.4)
abline(h=seq(0,1,by=0.1),lty=3,col=gray(0.5))
abline(v=seq(0,1400,by=200),lty=3,col=gray(0.5))
plot(sort(log10(Ba)),((1:n)-0.5)/n,pch=3,cex=0.8,
main="",xlab="Ba [mg/kg]",ylab="Probability",cex.lab=1,xaxt="n",cex.lab=1.4)
axis(1,at=log10(alog<-sort(c((10^(-50:50))%*%t(c(5,10))))),labels=alog)
abline(h=seq(0,1,by=0.1),lty=3,col=gray(0.5))
abline(v=log10(alog),lty=3,col=gray(0.5))


###################################################
### code chunk number 2: fig1
###################################################
library(StatDA)
data(chorizon)
Ba=chorizon[,"Ba"]
n=length(Ba)
par(mfcol=c(2,2),mar=c(4,4,2,2))
edaplotlog(Ba,H.freq=F,box=T,H.breaks=30,S.pch=3,S.cex=0.5,D.lwd=1.5,P.log=F,
  P.main="",P.xlab="Ba [mg/kg]",P.ylab="Density",B.pch=3,B.cex=0.5,B.log=TRUE)
edaplot(log10(Ba),H.freq=F,box=T,S.pch=3,S.cex=0.5,D.lwd=1.5,P.ylab="Density",
   P.log=T,P.logfine=c(5,10),P.main="",P.xlab="Ba [mg/kg]",B.pch=3,B.cex=0.5)
plot(sort(Ba),((1:n)-0.5)/n,pch=3,cex=0.8,
main="",xlab="Ba [mg/kg]",ylab="Probability",cex.lab=1,cex.lab=1.4)
abline(h=seq(0,1,by=0.1),lty=3,col=gray(0.5))
abline(v=seq(0,1400,by=200),lty=3,col=gray(0.5))
plot(sort(log10(Ba)),((1:n)-0.5)/n,pch=3,cex=0.8,
main="",xlab="Ba [mg/kg]",ylab="Probability",cex.lab=1,xaxt="n",cex.lab=1.4)
axis(1,at=log10(alog<-sort(c((10^(-50:50))%*%t(c(5,10))))),labels=alog)
abline(h=seq(0,1,by=0.1),lty=3,col=gray(0.5))
abline(v=log10(alog),lty=3,col=gray(0.5))


###################################################
### code chunk number 3: out
###################################################
"descriptive" <-
function (x)
{
    lab <- c("MIN", "Q1", "MEDIAN", "MEAN", "Q3", "MAX",
        "SD", "MAD", "CV %", "CVR %")
    #lab <- c("N", "Missings", "MIN", paste("Q_",quanlow <- c(0.02,0.05,0.1),sep=""),
    #    "Q1", "MEDIAN", "MEAN-log", "MEAN", "Q3", 
    #    paste("Q_",quanup <- c(0.9,0.95,0.98),sep=""), "MAX",
    #    "SD", "MAD", "IQR","CV", "CVR",
    #    "KS-norm", "SW-norm", "KS-lognorm", "SW-lognorm")
    if (missing(x)) {
        return(lab)
    }
    temp <- rep(0, length(lab))
    xt <- x[!is.na(x)]
    ix <- order(xt)
    n <- length(xt)
    if (!is.numeric(xt) || all(is.na(x))) {
        #return(c(n, rep(NA, length(lab) - 2), length(x) - length(xt)))
        return(c(rep(NA, length(lab) - 3), length(x) - length(xt)))
    }
    if (n == 1) {
        #return(c(n, xt[1], NA, rep(xt[1], 5), length(x) - length(xt)))
        return(c(NA, rep(xt[1], 5), length(x) - length(xt)))
    }
    else {
        #return(c(n, length(x) - length(xt), min(xt), quantile(xt,quanlow),
        return(c(min(xt), quantile(xt,0.25), median(xt), mean(xt), quantile(xt,0.75), max(xt),
            sd(xt), mad(xt), sd(xt)/mean(xt)*100, mad(xt)/median(xt)*100))
    }
}

"sumstats" <-
function (x, by)
{
    if (!missing(by)) {
        x <- cat.to.list(c(x), by)
    }
    if (!is.list(x) & !is.matrix(x))
        x <- matrix(x, ncol = 1)
    if (is.list(x)) {
        nrow <- length(x)
        out <- matrix(NA, nrow = nrow, ncol = length(descriptive()))
        dimnames(out) <- list(names(x), descriptive())
        for (j in (1:nrow)) {
            if (is.numeric(x[[j]])) {
                out[j,] <- descriptive(x[[j]])
            }
        }
        return(out)
    }
    if (is.matrix(x)) {
        nr <- ncol(x)
        out <- matrix(NA, nrow = nr, ncol = length(descriptive()))
        dimnames(out) <- list(dimnames(x)[[2]], descriptive())
        for (j in (1:nr)) {
            out[j, ] <- descriptive(x[, j])
        }
        return(out)
    }
}
data(moss)
out<-sumstats(moss[,20:58])



###################################################
### code chunk number 4: tab1
###################################################
xtable::xtable(out,caption="Summary statistics for selected elements of the Kola moss data.", label="tab1")




###################################################
### code chunk number 5: fig2plot
###################################################
par(mfrow=c(1,1))
data(kola.background)
el=chorizon[,"As"]
X=chorizon[,"XCOO"]
Y=chorizon[,"YCOO"]

plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)

bubbleFIN(X,Y,el,S=9,s=2,plottitle="",legendtitle="As [mg/kg]", text.cex=0.60,
	legtitle.cex=0.70,ndigits=2)

text(min(X)+diff(range(X))*5/7,max(Y),"Exponentially",cex=0.70)
text(min(X)+diff(range(X))*5/7,max(Y)-diff(range(Y))/25,"growing dots",cex=0.70)

scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=4e4,sizetext=0.8)

Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)


###################################################
### code chunk number 6: fig2
###################################################
par(mfrow=c(1,1))
data(kola.background)
el=chorizon[,"As"]
X=chorizon[,"XCOO"]
Y=chorizon[,"YCOO"]

plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)

bubbleFIN(X,Y,el,S=9,s=2,plottitle="",legendtitle="As [mg/kg]", text.cex=0.60,
	legtitle.cex=0.70,ndigits=2)

text(min(X)+diff(range(X))*5/7,max(Y),"Exponentially",cex=0.70)
text(min(X)+diff(range(X))*5/7,max(Y)-diff(range(Y))/25,"growing dots",cex=0.70)

scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=4e4,sizetext=0.8)

Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)


###################################################
### code chunk number 7: fig3plot
###################################################
el=chorizon[,"As"]
X=chorizon[,"XCOO"]
Y=chorizon[,"YCOO"]
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)

SymbLegend(X,Y,el,type="percentile",qutiles<-c(0,0.05,0.25,0.75,0.95,1),symbtype="EDA",symbmagn=0.8,
leg.position="topright",leg.title="As [mg/kg]",leg.title.cex=0.8,leg.round=2,leg.wid=4,leg.just="right")

text(min(X)+diff(range(X))*4/7,max(Y),paste(qutiles*100,collapse=","),cex=0.8)
text(min(X)+diff(range(X))*4/7,max(Y)-diff(range(Y))/25,"Percentiles",cex=0.8)

scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)


###################################################
### code chunk number 8: fig3
###################################################
el=chorizon[,"As"]
X=chorizon[,"XCOO"]
Y=chorizon[,"YCOO"]
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)

SymbLegend(X,Y,el,type="percentile",qutiles<-c(0,0.05,0.25,0.75,0.95,1),symbtype="EDA",symbmagn=0.8,
leg.position="topright",leg.title="As [mg/kg]",leg.title.cex=0.8,leg.round=2,leg.wid=4,leg.just="right")

text(min(X)+diff(range(X))*4/7,max(Y),paste(qutiles*100,collapse=","),cex=0.8)
text(min(X)+diff(range(X))*4/7,max(Y)-diff(range(Y))/25,"Percentiles",cex=0.8)

scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)


###################################################
### code chunk number 9: fig4plot
###################################################
X=chorizon[,"XCOO"]
Y=chorizon[,"YCOO"]
el=log10(chorizon[,"As"])
data(kola.background)
data(bordersKola)

plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")

SmoothLegend(X,Y,el,resol=200,type="contin",whichcol="gray",
    qutiles=c(0,0.05,0.25,0.50,0.75,0.90,0.95,1),borders="bordersKola",
    leg.xpos.min=7.8e5,leg.xpos.max=8.0e5,leg.ypos.min=77.6e5,leg.ypos.max=78.7e5,
    leg.title="mg/kg", leg.title.cex=0.7, leg.numb.cex=0.7, leg.round=2,leg.wid=4,
    leg.numb.xshift=0.7e5,leg.perc.xshift=0.4e5,leg.perc.yshift=0.2e5,tit.xshift=0.35e5)

plotbg(map.col=c("gray","gray","gray","gray"),map.lwd=c(1,1,1,1),add.plot=T)

text(min(X)+diff(range(X))*4/7,max(Y),"As",cex=1)
text(min(X)+diff(range(X))*4/7,max(Y)-diff(range(Y))/28,"in C-horizon",cex=0.8)

scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)


###################################################
### code chunk number 10: fig4
###################################################
X=chorizon[,"XCOO"]
Y=chorizon[,"YCOO"]
el=log10(chorizon[,"As"])
data(kola.background)
data(bordersKola)

plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")

SmoothLegend(X,Y,el,resol=200,type="contin",whichcol="gray",
    qutiles=c(0,0.05,0.25,0.50,0.75,0.90,0.95,1),borders="bordersKola",
    leg.xpos.min=7.8e5,leg.xpos.max=8.0e5,leg.ypos.min=77.6e5,leg.ypos.max=78.7e5,
    leg.title="mg/kg", leg.title.cex=0.7, leg.numb.cex=0.7, leg.round=2,leg.wid=4,
    leg.numb.xshift=0.7e5,leg.perc.xshift=0.4e5,leg.perc.yshift=0.2e5,tit.xshift=0.35e5)

plotbg(map.col=c("gray","gray","gray","gray"),map.lwd=c(1,1,1,1),add.plot=T)

text(min(X)+diff(range(X))*4/7,max(Y),"As",cex=1)
text(min(X)+diff(range(X))*4/7,max(Y)-diff(range(Y))/28,"in C-horizon",cex=0.8)

scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)


###################################################
### code chunk number 11: fig5imp
###################################################
X=chorizon[,"XCOO"]/1000
Y=chorizon[,"YCOO"]/1000
el=chorizon[,"As"]
vario.b <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300)
vario.c <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300, op="cloud")
vario.bc <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300, bin.cloud=TRUE)
vario.s <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300, op="sm", band=10)
vario.0 <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300, dir=0, tol=pi/8)
vario.90 <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300, dir=pi/4, tol=pi/8)
vario.45 <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300, dir=pi/8, tol=pi/8)
vario.120 <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300, dir=3*pi/8, tol=pi/8)
data(res.eyefit.As_C)
v5 <- variofit(vario.b,res.eyefit.As_C,cov.model="spherical",max.dist=300)


###################################################
### code chunk number 12: fig5plot
###################################################
par(mfrow=c(1,2),mar=c(4,4,2,2))
plot(0,0,xlab="Distance [km]",ylab="Semivariogram",xlim=c(0,300),ylim=c(0,1.6),
     cex.lab=1.2,type="n")
title("As in C-horizon",cex.main=1.2)
lines(vario.b,pch=1,type="p")
lines(vario.0,lty=1,pch=2)
lines(vario.90,lty=2,pch=3)
lines(vario.45,lty=3,pch=4)
lines(vario.120,lty=4,pch=5)
legend("bottomright", legend=c("omnidirectional","N-S","E-W","NW-SE","NE-SW"), pch=c(1,2,3,4,5))
data(res.eyefit.As_C)
plot(0,0,xlab="Distance [km]",ylab="Semivariogram",xlim=c(0,300),ylim=c(0,1.6),
     cex.lab=1.2,type="n")
title("As in C-horizon",cex.main=1.2)
lines(vario.b,pch=1,type="p")
lines(v5,col=1,lwd=2)
r=v5$cov.pars[2]
n=v5$nugget
s=v5$cov.pars[1]
arrows(0,0,0,n,length=0.08,code=3,col=gray(0.6))
text(2,n/2,paste("Nugget variance =",round(n,2)),cex=0.9,pos=4)
abline(h=n,col=gray(0.6),lty=2)
arrows(300,n,300,n+s,length=0.08,code=3,col=gray(0.6))
text(298,n+s/2,paste("Sill =",round(s,2)),cex=0.9,pos=2)
arrows(0,1.3,r,1.3,length=0.08,code=3,col=gray(0.6))
text(r/2,1.3,paste("Range =",round(r,0)),cex=0.9,pos=3)


###################################################
### code chunk number 13: fig5
###################################################
X=chorizon[,"XCOO"]/1000
Y=chorizon[,"YCOO"]/1000
el=chorizon[,"As"]
vario.b <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300)
vario.c <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300, op="cloud")
vario.bc <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300, bin.cloud=TRUE)
vario.s <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300, op="sm", band=10)
vario.0 <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300, dir=0, tol=pi/8)
vario.90 <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300, dir=pi/4, tol=pi/8)
vario.45 <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300, dir=pi/8, tol=pi/8)
vario.120 <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300, dir=3*pi/8, tol=pi/8)
data(res.eyefit.As_C)
v5 <- variofit(vario.b,res.eyefit.As_C,cov.model="spherical",max.dist=300)
par(mfrow=c(1,2),mar=c(4,4,2,2))
plot(0,0,xlab="Distance [km]",ylab="Semivariogram",xlim=c(0,300),ylim=c(0,1.6),
     cex.lab=1.2,type="n")
title("As in C-horizon",cex.main=1.2)
lines(vario.b,pch=1,type="p")
lines(vario.0,lty=1,pch=2)
lines(vario.90,lty=2,pch=3)
lines(vario.45,lty=3,pch=4)
lines(vario.120,lty=4,pch=5)
legend("bottomright", legend=c("omnidirectional","N-S","E-W","NW-SE","NE-SW"), pch=c(1,2,3,4,5))
data(res.eyefit.As_C)
plot(0,0,xlab="Distance [km]",ylab="Semivariogram",xlim=c(0,300),ylim=c(0,1.6),
     cex.lab=1.2,type="n")
title("As in C-horizon",cex.main=1.2)
lines(vario.b,pch=1,type="p")
lines(v5,col=1,lwd=2)
r=v5$cov.pars[2]
n=v5$nugget
s=v5$cov.pars[1]
arrows(0,0,0,n,length=0.08,code=3,col=gray(0.6))
text(2,n/2,paste("Nugget variance =",round(n,2)),cex=0.9,pos=4)
abline(h=n,col=gray(0.6),lty=2)
arrows(300,n,300,n+s,length=0.08,code=3,col=gray(0.6))
text(298,n+s/2,paste("Sill =",round(s,2)),cex=0.9,pos=2)
arrows(0,1.3,r,1.3,length=0.08,code=3,col=gray(0.6))
text(r/2,1.3,paste("Range =",round(r,0)),cex=0.9,pos=3)



###################################################
### code chunk number 14: fig6plot
###################################################
par(mfrow=c(1,1))
X=chorizon[,"XCOO"]
Y=chorizon[,"YCOO"]
xwid=diff(range(X))/12e4
ywid=diff(range(Y))/12e4
el=chorizon[,"As"]
vario.b <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300000)
data(res.eyefit.As_C_m)
data(bordersKola)
v5 <- variofit(vario.b,res.eyefit.As_C_m,cov.model="spherical",max.dist=300000)
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
KrigeLegend(X,Y,el,resol=25,vario=v5,type="percentile",whichcol="gray",
    qutiles=c(0,0.05,0.25,0.50,0.75,0.90,0.95,1),borders="bordersKola",
  leg.xpos.min=7.8e5,leg.xpos.max=8.0e5,leg.ypos.min=77.6e5,leg.ypos.max=78.7e5,
    leg.title="mg/kg", leg.title.cex=0.7, leg.numb.cex=0.7, leg.round=2,
    leg.numb.xshift=0.7e5,leg.perc.xshift=0.4e5,leg.perc.yshift=0.2e5,tit.xshift=0.35e5)
plotbg(map.col=c("gray","gray","gray","gray"),map.lwd=c(1,1,1,1),add.plot=T)
text(min(X)+diff(range(X))*4/7,max(Y),"As",cex=1)
text(min(X)+diff(range(X))*4/7,max(Y)-diff(range(Y))/28,"in C-horizon",cex=0.8)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)


###################################################
### code chunk number 15: fig6
###################################################
par(mfrow=c(1,1))
X=chorizon[,"XCOO"]
Y=chorizon[,"YCOO"]
xwid=diff(range(X))/12e4
ywid=diff(range(Y))/12e4
el=chorizon[,"As"]
vario.b <- variog(coords=cbind(X,Y), data=el, lambda=0, max.dist=300000)
data(res.eyefit.As_C_m)
data(bordersKola)
v5 <- variofit(vario.b,res.eyefit.As_C_m,cov.model="spherical",max.dist=300000)
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
KrigeLegend(X,Y,el,resol=25,vario=v5,type="percentile",whichcol="gray",
    qutiles=c(0,0.05,0.25,0.50,0.75,0.90,0.95,1),borders="bordersKola",
  leg.xpos.min=7.8e5,leg.xpos.max=8.0e5,leg.ypos.min=77.6e5,leg.ypos.max=78.7e5,
    leg.title="mg/kg", leg.title.cex=0.7, leg.numb.cex=0.7, leg.round=2,
    leg.numb.xshift=0.7e5,leg.perc.xshift=0.4e5,leg.perc.yshift=0.2e5,tit.xshift=0.35e5)
plotbg(map.col=c("gray","gray","gray","gray"),map.lwd=c(1,1,1,1),add.plot=T)
text(min(X)+diff(range(X))*4/7,max(Y),"As",cex=1)
text(min(X)+diff(range(X))*4/7,max(Y)-diff(range(Y))/28,"in C-horizon",cex=0.8)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)


###################################################
### code chunk number 16: fig7plot
###################################################
data(moss)
par(mfrow=c(1,2),mar=c(4,4,2,2))

plotbg(map.col=c("gray","gray","gray","gray"), xlab="UTM east [m]", ylab="UTM north [m]",
cex.lab=1.2)
X=moss[,"XCOO"]
Y=moss[,"YCOO"]
points(X[Y<7600000 & Y>7500000],Y[Y<7600000 & Y>7500000],pch=3,cex=0.7)
x=(X[Y<7600000 & Y>7500000]-753970)/1000
y=log10(moss[Y<7600000 & Y>7500000,"Cu"])
plot(x,y,xlab="Distance from Monchegorsk [km]",ylab="Cu in Moss [mg/kg]", yaxt="n",
     pch=3,cex=0.7, cex.lab=1.2)
axis(2,at=log10(a<-sort(c((10^(-50:50))%*%t(c(2,5,10))))),labels=a)
lines(smooth.spline(x,y),col=1,lwd=1.3)


###################################################
### code chunk number 17: fig7
###################################################
data(moss)
par(mfrow=c(1,2),mar=c(4,4,2,2))

plotbg(map.col=c("gray","gray","gray","gray"), xlab="UTM east [m]", ylab="UTM north [m]",
cex.lab=1.2)
X=moss[,"XCOO"]
Y=moss[,"YCOO"]
points(X[Y<7600000 & Y>7500000],Y[Y<7600000 & Y>7500000],pch=3,cex=0.7)
x=(X[Y<7600000 & Y>7500000]-753970)/1000
y=log10(moss[Y<7600000 & Y>7500000,"Cu"])
plot(x,y,xlab="Distance from Monchegorsk [km]",ylab="Cu in Moss [mg/kg]", yaxt="n",
     pch=3,cex=0.7, cex.lab=1.2)
axis(2,at=log10(a<-sort(c((10^(-50:50))%*%t(c(2,5,10))))),labels=a)
lines(smooth.spline(x,y),col=1,lwd=1.3)


###################################################
### code chunk number 18: fig8plot
###################################################
X=moss[,"XCOO"]
Y=moss[,"YCOO"]
Cu=moss[,"Cu"]
loc = list(x=c(590565.1,511872.8,779702.8,861156.3),
           y=c(7824652,7769344,7428054,7510341))
loc.in=in.polygon(X,Y,loc$x,loc$y) # observations in polygone
ref=c(542245.3,7803068) # reference point
par(mfrow=c(1,2),mar=c(4,4,2,2))
plotbg(map.col=c("gray","gray","gray","gray"), xlab="UTM east [m]", ylab="UTM north [m]",
cex.lab=1.2)
points(X[!loc.in],Y[!loc.in],pch=16,cex=0.3)
points(X[loc.in],Y[loc.in],pch=3,cex=0.7)
points(ref[1],ref[2],pch=16,cex=1)
polygon(loc$x,loc$y,border=1,lwd=1.3)
distanc=sqrt((X[loc.in]-ref[1])^2+(Y[loc.in]-ref[2])^2)
plot(distanc/1000,log10(Cu[loc.in]),xlab="Distance from reference point [km]",
       ylab="Cu in Moss [mg/kg]", yaxt="n", pch=3,cex=0.7, cex.lab=1.2)
axis(2,at=log10(a<-sort(c((10^(-50:50))%*%t(c(2,5,10))))),labels=a)
lines(smooth.spline(distanc/1000,log10(Cu[loc.in]),spar=0.8),col=1,lwd=1.3)


###################################################
### code chunk number 19: fig8
###################################################
X=moss[,"XCOO"]
Y=moss[,"YCOO"]
Cu=moss[,"Cu"]
loc = list(x=c(590565.1,511872.8,779702.8,861156.3),
           y=c(7824652,7769344,7428054,7510341))
loc.in=in.polygon(X,Y,loc$x,loc$y) # observations in polygone
ref=c(542245.3,7803068) # reference point
par(mfrow=c(1,2),mar=c(4,4,2,2))
plotbg(map.col=c("gray","gray","gray","gray"), xlab="UTM east [m]", ylab="UTM north [m]",
cex.lab=1.2)
points(X[!loc.in],Y[!loc.in],pch=16,cex=0.3)
points(X[loc.in],Y[loc.in],pch=3,cex=0.7)
points(ref[1],ref[2],pch=16,cex=1)
polygon(loc$x,loc$y,border=1,lwd=1.3)
distanc=sqrt((X[loc.in]-ref[1])^2+(Y[loc.in]-ref[2])^2)
plot(distanc/1000,log10(Cu[loc.in]),xlab="Distance from reference point [km]",
       ylab="Cu in Moss [mg/kg]", yaxt="n", pch=3,cex=0.7, cex.lab=1.2)
axis(2,at=log10(a<-sort(c((10^(-50:50))%*%t(c(2,5,10))))),labels=a)
lines(smooth.spline(distanc/1000,log10(Cu[loc.in]),spar=0.8),col=1,lwd=1.3)


###################################################
### code chunk number 20: fig9plot
###################################################
par(mfrow=c(1,1))
x=moss[,c("Al","Fe","Mn")]
ternary(x,grid=TRUE,pch=3,cex=0.7,col=1)
text(0.21,0.26,"563")
text(0.1,0.05,"726")
text(0.21,0.06,"325")


###################################################
### code chunk number 21: fig9
###################################################
par(mfrow=c(1,1))
x=moss[,c("Al","Fe","Mn")]
ternary(x,grid=TRUE,pch=3,cex=0.7,col=1)
text(0.21,0.26,"563")
text(0.1,0.05,"726")
text(0.21,0.06,"325")


###################################################
### code chunk number 22: fig10plot
###################################################
data(ohorizon)
data(bhorizon)

nam="Al"
M=density(log10(moss[,nam]))
O=density(log10(ohorizon[,nam]))
B=density(log10(bhorizon[,nam]))
C=density(log10(chorizon[,nam]))
plot(0,0,xlim=c(min(M$x,O$x,B$x,C$x),max(M$x,O$x,B$x,C$x)),
  ylim=c(min(M$y,O$y,B$y,C$y),max(M$y,O$y,B$y,C$y)),
  main="",xlab=paste(nam," [mg/kg]",sep=""),ylab="Density",
  xaxt="n",cex.lab=1.2,type="n")
axis(1,at=log10(a<-sort(c((10^(-50:50))%*%t(c(2,5,10))))),labels=a)
lines(M,lty=1)
lines(O,lty=2)
lines(B,lty=4)
lines(C,lty=5)
text(1.5,1.5,"Moss",pos=4)
text(2.5,1.45,"O-horizon",pos=4)
text(4.6,0.9,"B-horizon",pos=4,srt=90)
text(3.7,1.4,"C-horizon",pos=4,srt=90)


###################################################
### code chunk number 23: fig10
###################################################
data(ohorizon)
data(bhorizon)

nam="Al"
M=density(log10(moss[,nam]))
O=density(log10(ohorizon[,nam]))
B=density(log10(bhorizon[,nam]))
C=density(log10(chorizon[,nam]))
plot(0,0,xlim=c(min(M$x,O$x,B$x,C$x),max(M$x,O$x,B$x,C$x)),
  ylim=c(min(M$y,O$y,B$y,C$y),max(M$y,O$y,B$y,C$y)),
  main="",xlab=paste(nam," [mg/kg]",sep=""),ylab="Density",
  xaxt="n",cex.lab=1.2,type="n")
axis(1,at=log10(a<-sort(c((10^(-50:50))%*%t(c(2,5,10))))),labels=a)
lines(M,lty=1)
lines(O,lty=2)
lines(B,lty=4)
lines(C,lty=5)
text(1.5,1.5,"Moss",pos=4)
text(2.5,1.45,"O-horizon",pos=4)
text(4.6,0.9,"B-horizon",pos=4,srt=90)
text(3.7,1.4,"C-horizon",pos=4,srt=90)


###################################################
### code chunk number 24: fig11plot
###################################################
nam="Al"
M=log10(moss[,nam])
O=log10(ohorizon[,nam])
B=log10(bhorizon[,nam])
C=log10(chorizon[,nam])
qpplot.das(M,qdist=qnorm,xlab=paste(nam," [mg/kg]",sep=""),
ylab="Probability [%]", pch=3,cex=0.7, logx=TRUE,
logfinetick=c(2,5,10),logfinelab=c(2,5,10),line=F,cex.lab=1.2,
xlim=c(min(M,O,B,C),max(M,O,B,C)))
points(sort(O),qnorm(ppoints(length(O))),pch=4,cex=0.7)
points(sort(B),qnorm(ppoints(length(B))),pch=1,cex=0.7)
points(sort(C),qnorm(ppoints(length(C))),pch=22,cex=0.7)
legend("topleft",legend=c("Moss","O-horizon","B-horizon","C-horizon"),
       pch=c(3,4,1,22),bg="white")


###################################################
### code chunk number 25: fig11
###################################################
nam="Al"
M=log10(moss[,nam])
O=log10(ohorizon[,nam])
B=log10(bhorizon[,nam])
C=log10(chorizon[,nam])
qpplot.das(M,qdist=qnorm,xlab=paste(nam," [mg/kg]",sep=""),
ylab="Probability [%]", pch=3,cex=0.7, logx=TRUE,
logfinetick=c(2,5,10),logfinelab=c(2,5,10),line=F,cex.lab=1.2,
xlim=c(min(M,O,B,C),max(M,O,B,C)))
points(sort(O),qnorm(ppoints(length(O))),pch=4,cex=0.7)
points(sort(B),qnorm(ppoints(length(B))),pch=1,cex=0.7)
points(sort(C),qnorm(ppoints(length(C))),pch=22,cex=0.7)
legend("topleft",legend=c("Moss","O-horizon","B-horizon","C-horizon"),
       pch=c(3,4,1,22),bg="white")


###################################################
### code chunk number 26: fig12plot
###################################################
nam="Al"
M=moss[,nam]
O=ohorizon[,nam]
B=bhorizon[,nam]
C=chorizon[,nam]
boxplotlog(C,B,O,M,horizontal=T,xlab=paste(nam," [mg/kg]",sep=""),cex=0.7,cex.axis=1.2,
   log="x",cex.lab=1.2,pch=3,names=c("C-hor","B-hor","O-hor","Moss"),xaxt="n")
axis(1,at=(a<-sort(c((10^(-50:50))%*%t(c(2,5,10))))),labels=a)


###################################################
### code chunk number 27: fig12
###################################################
nam="Al"
M=moss[,nam]
O=ohorizon[,nam]
B=bhorizon[,nam]
C=chorizon[,nam]
boxplotlog(C,B,O,M,horizontal=T,xlab=paste(nam," [mg/kg]",sep=""),cex=0.7,cex.axis=1.2,
   log="x",cex.lab=1.2,pch=3,names=c("C-hor","B-hor","O-hor","Moss"),xaxt="n")
axis(1,at=(a<-sort(c((10^(-50:50))%*%t(c(2,5,10))))),labels=a)


###################################################
### code chunk number 28: fig13plot
###################################################
x=chorizon[,c("Ca","Cu","Mg","Na","P","Sr","Zn")]
par(mfrow=c(1,1),mar=c(4,4,2,0))
R=robustbase::covMcd(log10(x),cor=T)$cor
P=cor(log10(x))
CorCompare(R,P,labels1=dimnames(x)[[2]],labels2=dimnames(x)[[2]],
method1="Robust",method2="Pearson",ndigits=2, cex.label=1.2)


###################################################
### code chunk number 29: fig13
###################################################
x=chorizon[,c("Ca","Cu","Mg","Na","P","Sr","Zn")]
par(mfrow=c(1,1),mar=c(4,4,2,0))
R=robustbase::covMcd(log10(x),cor=T)$cor
P=cor(log10(x))
CorCompare(R,P,labels1=dimnames(x)[[2]],labels2=dimnames(x)[[2]],
method1="Robust",method2="Pearson",ndigits=2, cex.label=1.2)


###################################################
### code chunk number 30: fig14plot
###################################################
X=ohorizon[,"XCOO"]
Y=ohorizon[,"YCOO"]
el=log10(ohorizon[,c("Co","Cu","Ni","Rb","Bi","Na","Sr")])
data(kola.background)
sel <- c(3,8,22, 29, 32, 35, 43, 69, 73 ,93,109,129,130,134,168,181,183,205,211,
      218,237,242,276,292,297,298,345,346,352,372,373,386,408,419,427,441,446,490,
      516,535,551,556,558,564,577,584,601,612,617)
x=el[sel,]
dimnames(x)[[1]]=ohorizon[sel,1]
xwid=diff(range(X))/12e4
ywid=diff(range(Y))/12e4
par(mfrow=c(1,2),mar=c(1.5,1.5,1.5,1.5))
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n",
   xlim=c(360000,max(X)))
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)
tree(x,locations=cbind(X[sel],Y[sel]),len=700,key.loc=c(793000,7760000),leglen=1500,
     cex=0.75, add=T, leglh=6,lh=30,wmax=120,wmin=30,labels=NULL)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)
par(mar=c(.1,.1,.1,.5))
tree(x,key.loc=c(15,0),len=0.022, lh=30,leglh=4,
    wmax=120,wmin=30, leglen=0.05, ncol=8, cex=0.75)


###################################################
### code chunk number 31: fig14
###################################################
X=ohorizon[,"XCOO"]
Y=ohorizon[,"YCOO"]
el=log10(ohorizon[,c("Co","Cu","Ni","Rb","Bi","Na","Sr")])
data(kola.background)
sel <- c(3,8,22, 29, 32, 35, 43, 69, 73 ,93,109,129,130,134,168,181,183,205,211,
      218,237,242,276,292,297,298,345,346,352,372,373,386,408,419,427,441,446,490,
      516,535,551,556,558,564,577,584,601,612,617)
x=el[sel,]
dimnames(x)[[1]]=ohorizon[sel,1]
xwid=diff(range(X))/12e4
ywid=diff(range(Y))/12e4
par(mfrow=c(1,2),mar=c(1.5,1.5,1.5,1.5))
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n",
   xlim=c(360000,max(X)))
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)
tree(x,locations=cbind(X[sel],Y[sel]),len=700,key.loc=c(793000,7760000),leglen=1500,
     cex=0.75, add=T, leglh=6,lh=30,wmax=120,wmin=30,labels=NULL)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)
par(mar=c(.1,.1,.1,.5))
tree(x,key.loc=c(15,0),len=0.022, lh=30,leglh=4,
    wmax=120,wmin=30, leglen=0.05, ncol=8, cex=0.75)


###################################################
### code chunk number 32: fig15plot
###################################################
X=moss[,"XCOO"]
Y=moss[,"YCOO"]
xwid=diff(range(X))/12e4
ywid=diff(range(Y))/12e4
par(mfrow=c(1,2),mar=c(1.5,1.5,1.5,1.5))
el=c("Ag","As","Bi","Cd","Co","Cu","Ni")
dat=log10(moss[,el])
res <- plotmvoutlier(cbind(X,Y),dat,symb=F,map.col=c("grey","grey","grey","grey"),
       map.lwd=c(1,1,1,1),
       xlab="",ylab="",frame.plot=FALSE,xaxt="n",yaxt="n")    
legend("topright",pch=c(3,21),pt.cex=c(0.7,0.2), legend=c("outliers",
"non-outliers"), cex=0.8)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)
res <- plotmvoutlier(cbind(X,Y),dat,symb=T,bw=T,map.col=c("grey","grey","grey","grey"),
       map.lwd=c(1,1,1,1), lcex.fac=0.6,
       xlab="",ylab="",frame.plot=FALSE,xaxt="n",yaxt="n")    
perc.ypos=seq(from=77.4e5,to=78.7e5,length=6)
sym.ypos=perc.ypos[-6]+diff(perc.ypos)/2
points(rep(8.2e5,5),sym.ypos,pch=c(1, 1, 16, 3, 3),cex=c(1.5, 1, 0.5, 1, 1.5)*0.6)
text(rep(8.0e5,6),perc.ypos[-5],round(100*c(0,0.25,0.5,0.75,1),2),pos=2,cex=0.7)
text(8.0e5,perc.ypos[5],"break",cex=0.72,pos=2)
text(7.5e5,79e5,"Percentile",cex=0.8)
text(8.4e5,79e5,"Symbol",cex=0.8)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)


###################################################
### code chunk number 33: fig15
###################################################
X=moss[,"XCOO"]
Y=moss[,"YCOO"]
xwid=diff(range(X))/12e4
ywid=diff(range(Y))/12e4
par(mfrow=c(1,2),mar=c(1.5,1.5,1.5,1.5))
el=c("Ag","As","Bi","Cd","Co","Cu","Ni")
dat=log10(moss[,el])
res <- plotmvoutlier(cbind(X,Y),dat,symb=F,map.col=c("grey","grey","grey","grey"),
       map.lwd=c(1,1,1,1),
       xlab="",ylab="",frame.plot=FALSE,xaxt="n",yaxt="n")    
legend("topright",pch=c(3,21),pt.cex=c(0.7,0.2), legend=c("outliers",
"non-outliers"), cex=0.8)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)
res <- plotmvoutlier(cbind(X,Y),dat,symb=T,bw=T,map.col=c("grey","grey","grey","grey"),
       map.lwd=c(1,1,1,1), lcex.fac=0.6,
       xlab="",ylab="",frame.plot=FALSE,xaxt="n",yaxt="n")    
perc.ypos=seq(from=77.4e5,to=78.7e5,length=6)
sym.ypos=perc.ypos[-6]+diff(perc.ypos)/2
points(rep(8.2e5,5),sym.ypos,pch=c(1, 1, 16, 3, 3),cex=c(1.5, 1, 0.5, 1, 1.5)*0.6)
text(rep(8.0e5,6),perc.ypos[-5],round(100*c(0,0.25,0.5,0.75,1),2),pos=2,cex=0.7)
text(8.0e5,perc.ypos[5],"break",cex=0.72,pos=2)
text(7.5e5,79e5,"Percentile",cex=0.8)
text(8.4e5,79e5,"Symbol",cex=0.8)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)


###################################################
### code chunk number 34: fig16plot
###################################################
par(mfrow=c(1,1))
X=moss[,"XCOO"]
Y=moss[,"YCOO"]
var=c("Ag","Al","As","B","Ba","Bi","Ca","Cd","Co","Cr","Cu","Fe","Hg","K","Mg","Mn","Mo",
      "Na","Ni","P","Pb","Rb","S","Sb","Si","Sr","Th","Tl","U","V","Zn")
x=moss[,var]
ilr <- function(x){
x.ilr=matrix(NA,nrow=nrow(x),ncol=ncol(x)-1)
  for (i in 1:ncol(x.ilr)){
    x.ilr[,i]=sqrt((i)/(i+1))*log(((apply(as.matrix(x[,1:i]), 1, prod))^(1/i))/(x[,i+1]))
  }
return(x.ilr)
}
invilr <- function(x.ilr,x.ori){
y=matrix(0,nrow=nrow(x.ilr),ncol=ncol(x.ilr)+1)
for (i in 1:(ncol(y)-1)){
    for (j in i:ncol(x.ilr)){
      y[,i]=y[,i]+x.ilr[,j]/sqrt(j*(j+1))
    }
}
for (i in 2:ncol(y)){
   y[,i]=y[,i]-sqrt((i-1)/i)*x.ilr[,i-1]
}
yexp=exp(y)
x.back=yexp/apply(yexp,1,sum)
return(x.back)
}
x2=ilr(x)
res2=princomp(x2,cor=TRUE)
x2.mcd=robustbase::covMcd(x2,cor=T)
res2rob=princomp(x2,covmat=x2.mcd,cor=T)
V=matrix(0,nrow=31,ncol=30)
for (i in 1:ncol(V)){
  V[1:i,i] <- 1/i
  V[i+1,i] <- (-1)
  V[,i] <- V[,i]*sqrt(i/(i+1))
}
xgeom=10^apply(log10(x),1,mean)
xlc1=x/xgeom
xlc=log10(xlc1)
reslc=princomp(xlc,cor=TRUE)
res2robback=res2rob
res2robback$loadings <- V%*%res2rob$loadings
res2robback$scores <- res2rob$scores%*%t(V)
dimnames(res2robback$loadings)[[1]] <- names(x)
res2robback$sdev <- apply(res2robback$scores,2,mad)
biplot(res2robback,xlab="PC1 (robust)",ylab="PC2 (robust)",col=c(gray(0.6),1),
       xlabs=rep("+",nrow(x)),cex=0.8)


###################################################
### code chunk number 35: fig16
###################################################
par(mfrow=c(1,1))
X=moss[,"XCOO"]
Y=moss[,"YCOO"]
var=c("Ag","Al","As","B","Ba","Bi","Ca","Cd","Co","Cr","Cu","Fe","Hg","K","Mg","Mn","Mo",
      "Na","Ni","P","Pb","Rb","S","Sb","Si","Sr","Th","Tl","U","V","Zn")
x=moss[,var]
ilr <- function(x){
x.ilr=matrix(NA,nrow=nrow(x),ncol=ncol(x)-1)
  for (i in 1:ncol(x.ilr)){
    x.ilr[,i]=sqrt((i)/(i+1))*log(((apply(as.matrix(x[,1:i]), 1, prod))^(1/i))/(x[,i+1]))
  }
return(x.ilr)
}
invilr <- function(x.ilr,x.ori){
y=matrix(0,nrow=nrow(x.ilr),ncol=ncol(x.ilr)+1)
for (i in 1:(ncol(y)-1)){
    for (j in i:ncol(x.ilr)){
      y[,i]=y[,i]+x.ilr[,j]/sqrt(j*(j+1))
    }
}
for (i in 2:ncol(y)){
   y[,i]=y[,i]-sqrt((i-1)/i)*x.ilr[,i-1]
}
yexp=exp(y)
x.back=yexp/apply(yexp,1,sum)
return(x.back)
}
x2=ilr(x)
res2=princomp(x2,cor=TRUE)
x2.mcd=robustbase::covMcd(x2,cor=T)
res2rob=princomp(x2,covmat=x2.mcd,cor=T)
V=matrix(0,nrow=31,ncol=30)
for (i in 1:ncol(V)){
  V[1:i,i] <- 1/i
  V[i+1,i] <- (-1)
  V[,i] <- V[,i]*sqrt(i/(i+1))
}
xgeom=10^apply(log10(x),1,mean)
xlc1=x/xgeom
xlc=log10(xlc1)
reslc=princomp(xlc,cor=TRUE)
res2robback=res2rob
res2robback$loadings <- V%*%res2rob$loadings
res2robback$scores <- res2rob$scores%*%t(V)
dimnames(res2robback$loadings)[[1]] <- names(x)
res2robback$sdev <- apply(res2robback$scores,2,mad)
biplot(res2robback,xlab="PC1 (robust)",ylab="PC2 (robust)",col=c(gray(0.6),1),
       xlabs=rep("+",nrow(x)),cex=0.8)


###################################################
### code chunk number 36: fig17plot
###################################################
el=c("Cu","Ni","Mg","Sr")
x=scale(log10(moss[,el]))
set.seed(100)
res=e1071::cmeans(x,4)
X=moss[,"XCOO"]
Y=moss[,"YCOO"]
xwid=diff(range(X))/12e4
ywid=diff(range(Y))/12e4
par(mfrow=c(2,2),mar=c(1.5,1.5,1.5,1.5))
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)
points(moss[,"XCOO"],moss[,"YCOO"],col=gray(1-res$mem[,1]),pch=15,cex=1)
text(752000,7880000,"Cluster 1",cex=1.1)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)
points(moss[,"XCOO"],moss[,"YCOO"],col=gray(1-res$mem[,2]),pch=15,cex=1)
text(752000,7880000,"Cluster 2",cex=1.1)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)
points(moss[,"XCOO"],moss[,"YCOO"],col=gray(1-res$mem[,3]),pch=15,cex=1)
text(752000,7880000,"Cluster 3",cex=1.1)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)
points(moss[,"XCOO"],moss[,"YCOO"],col=gray(1-res$mem[,4]),pch=15,cex=1)
text(752000,7880000,"Cluster 4",cex=1.1)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)


###################################################
### code chunk number 37: fig17
###################################################
el=c("Cu","Ni","Mg","Sr")
x=scale(log10(moss[,el]))
set.seed(100)
res=e1071::cmeans(x,4)
X=moss[,"XCOO"]
Y=moss[,"YCOO"]
xwid=diff(range(X))/12e4
ywid=diff(range(Y))/12e4
par(mfrow=c(2,2),mar=c(1.5,1.5,1.5,1.5))
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)
points(moss[,"XCOO"],moss[,"YCOO"],col=gray(1-res$mem[,1]),pch=15,cex=1)
text(752000,7880000,"Cluster 1",cex=1.1)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)
points(moss[,"XCOO"],moss[,"YCOO"],col=gray(1-res$mem[,2]),pch=15,cex=1)
text(752000,7880000,"Cluster 2",cex=1.1)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)
points(moss[,"XCOO"],moss[,"YCOO"],col=gray(1-res$mem[,3]),pch=15,cex=1)
text(752000,7880000,"Cluster 3",cex=1.1)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)
points(moss[,"XCOO"],moss[,"YCOO"],col=gray(1-res$mem[,4]),pch=15,cex=1)
text(752000,7880000,"Cluster 4",cex=1.1)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)


###################################################
### code chunk number 38: fig18plot
###################################################
attach(chorizon)
X=chorizon[,"XCOO"]
Y=chorizon[,"YCOO"]
xdat=log10(chorizon[,101:109]/chorizon[,110])
set.seed(101)
xdat.mcd=robustbase::covMcd(xdat,cor=T)
md=sqrt(xdat.mcd$mah)
crit=sqrt(qchisq(0.975,ncol(xdat)))
set.seed(104)
res=robustbase::ltsReg(log10(Be) ~ Al_XRF+Ca_XRF+Fe_XRF+K_XRF+Mg_XRF+Mn_XRF+Na_XRF+P_XRF+Si_XRF,
   data=xdat)
resid=res$residuals/res$scale
xwid=diff(range(X))/12e4
ywid=diff(range(Y))/12e4
par(mfrow=c(2,2),mar=c(4,4,2,2))
psymb=res$lts.wt
psymb[res$lts.wt==1] <- 1
psymb[res$lts.wt==0] <- 3
pcex=res$lts.wt
pcex[res$lts.wt==1] <- 1.3
pcex[res$lts.wt==0] <- 0.8
qqnorm(resid,xlab="Quantiles of standard normal distribution",ylab="Standardised LTS residuals",
   pch=psymb,cex=pcex,main="",cex.lab=1.2)
qqline(resid)
plot(res$fitted,resid,cex.lab=1.2,xlab="Fitted values",ylab="Standardised LTS residuals",type="n")
points(res$fitted[res$lts.wt==0],resid[res$lts.wt==0],cex=0.8,pch=3)
points(res$fitted[res$lts.wt==1],resid[res$lts.wt==1],cex=0.8,pch=1)
abline(h=0,col="grey",lty=2)
abline(h=c(-2.5,2.5),lty=3,cex=1.1)
symb.nor=16 #1
symb.resl=1
symb.resh=22
symb.goodl=3 #16
symb.badll=1 #16
symb.badlh=15
plot(md,resid,cex=0.5,pch=3,type="n",xlab="Robust Mahalanobis distances",
        ylab="Standardised LTS residuals", cex.lab=1.2)
abline(h=c(2.5,-2.5))
abline(v=crit)
md.resid=as.data.frame(cbind(md,resid))
points(md.resid[md<crit & abs(resid)<2.5,], cex=0.3,pch=symb.nor)
points(md.resid[md<crit & resid>=2.5,], cex=0.9,pch=symb.resh)
points(md.resid[md<crit & resid<=(-2.5),], cex=0.9,pch=symb.resl)
points(md.resid[md>=crit & abs(resid)<2.5,], cex=0.6,pch=symb.goodl)
points(md.resid[md>=crit & resid>=2.5,], cex=0.9,pch=symb.badlh)
points(md.resid[md>=crit & resid<=(-2.5),], cex=1.1,pch=symb.badll)
par(mar=c(1.5,1.5,1.5,1.5))
XY=cbind(X,Y)
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)
points(XY[md<crit & abs(resid)<2.5,], cex=0.3,pch=symb.nor)
points(XY[md<crit & resid>=2.5,], cex=0.9,pch=symb.resh)
points(XY[md<crit & resid<=(-2.5),], cex=0.9,pch=symb.resl)
points(XY[md>=crit & abs(resid)<2.5,], cex=0.6,pch=symb.goodl)
points(XY[md>=crit & resid>=2.5,], cex=0.9,pch=symb.badlh)
points(XY[md>=crit & resid<=(-2.5),], cex=1.1,pch=symb.badll)
legend("topright",pch=c(symb.nor,symb.resh,symb.resl,symb.goodl,symb.badlh,symb.badll),
  pt.cex=c(0.3,0.9,0.9,0.6,0.9,1.1), 
legend=c("Normal observations","High vertical outliers","Low vertical outliers",
"Good leverage points","High bad leverage points","Low bad leverage points"), cex=0.7)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)


###################################################
### code chunk number 39: fig18
###################################################
attach(chorizon)
X=chorizon[,"XCOO"]
Y=chorizon[,"YCOO"]
xdat=log10(chorizon[,101:109]/chorizon[,110])
set.seed(101)
xdat.mcd=robustbase::covMcd(xdat,cor=T)
md=sqrt(xdat.mcd$mah)
crit=sqrt(qchisq(0.975,ncol(xdat)))
set.seed(104)
res=robustbase::ltsReg(log10(Be) ~ Al_XRF+Ca_XRF+Fe_XRF+K_XRF+Mg_XRF+Mn_XRF+Na_XRF+P_XRF+Si_XRF,
   data=xdat)
resid=res$residuals/res$scale
xwid=diff(range(X))/12e4
ywid=diff(range(Y))/12e4
par(mfrow=c(2,2),mar=c(4,4,2,2))
psymb=res$lts.wt
psymb[res$lts.wt==1] <- 1
psymb[res$lts.wt==0] <- 3
pcex=res$lts.wt
pcex[res$lts.wt==1] <- 1.3
pcex[res$lts.wt==0] <- 0.8
qqnorm(resid,xlab="Quantiles of standard normal distribution",ylab="Standardised LTS residuals",
   pch=psymb,cex=pcex,main="",cex.lab=1.2)
qqline(resid)
plot(res$fitted,resid,cex.lab=1.2,xlab="Fitted values",ylab="Standardised LTS residuals",type="n")
points(res$fitted[res$lts.wt==0],resid[res$lts.wt==0],cex=0.8,pch=3)
points(res$fitted[res$lts.wt==1],resid[res$lts.wt==1],cex=0.8,pch=1)
abline(h=0,col="grey",lty=2)
abline(h=c(-2.5,2.5),lty=3,cex=1.1)
symb.nor=16 #1
symb.resl=1
symb.resh=22
symb.goodl=3 #16
symb.badll=1 #16
symb.badlh=15
plot(md,resid,cex=0.5,pch=3,type="n",xlab="Robust Mahalanobis distances",
        ylab="Standardised LTS residuals", cex.lab=1.2)
abline(h=c(2.5,-2.5))
abline(v=crit)
md.resid=as.data.frame(cbind(md,resid))
points(md.resid[md<crit & abs(resid)<2.5,], cex=0.3,pch=symb.nor)
points(md.resid[md<crit & resid>=2.5,], cex=0.9,pch=symb.resh)
points(md.resid[md<crit & resid<=(-2.5),], cex=0.9,pch=symb.resl)
points(md.resid[md>=crit & abs(resid)<2.5,], cex=0.6,pch=symb.goodl)
points(md.resid[md>=crit & resid>=2.5,], cex=0.9,pch=symb.badlh)
points(md.resid[md>=crit & resid<=(-2.5),], cex=1.1,pch=symb.badll)
par(mar=c(1.5,1.5,1.5,1.5))
XY=cbind(X,Y)
plot(X,Y,frame.plot=FALSE,xaxt="n",yaxt="n",xlab="",ylab="",type="n")
plotbg(map.col=c("gray","gray","gray","gray"),add.plot=T)
points(XY[md<crit & abs(resid)<2.5,], cex=0.3,pch=symb.nor)
points(XY[md<crit & resid>=2.5,], cex=0.9,pch=symb.resh)
points(XY[md<crit & resid<=(-2.5),], cex=0.9,pch=symb.resl)
points(XY[md>=crit & abs(resid)<2.5,], cex=0.6,pch=symb.goodl)
points(XY[md>=crit & resid>=2.5,], cex=0.9,pch=symb.badlh)
points(XY[md>=crit & resid<=(-2.5),], cex=1.1,pch=symb.badll)
legend("topright",pch=c(symb.nor,symb.resh,symb.resl,symb.goodl,symb.badlh,symb.badll),
  pt.cex=c(0.3,0.9,0.9,0.6,0.9,1.1), 
legend=c("Normal observations","High vertical outliers","Low vertical outliers",
"Good leverage points","High bad leverage points","Low bad leverage points"), cex=0.7)
scalebar(761309,7373050,861309,7363050,shifttext=-0.5,shiftkm=37e3,sizetext=0.8)
Northarrow(362602,7818750,362602,7878750,362602,7838750,Alength=0.15,Aangle=15,Alwd=1.3,Tcex=1.6)


