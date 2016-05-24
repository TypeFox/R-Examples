aBombHome="~/data/abomb"
load(file.path(aBombHome,"Hema87.RData"));
names(d)
d=d[d$py>0,] #remove two recs with zero py
d=d[d$kerma==1,] # take only kerma < 4 Gy
d$py=10^4*d$py

svc<-with(d,cut(sv,c(-1,.02,1,10)))
PYT<-with(d,tapply(py,list(svc,agexg),sum,na.rm=TRUE))
AGET<-with(d,tapply(agex,list(svc,agexg),mean,na.rm=TRUE))
CMLT<-with(d,tapply(CML,list(svc,agexg),sum,na.rm=TRUE)) 
(CMLt=cbind(apply(CMLT[,1:4],1,sum),apply(CMLT[,5:8],1,sum),apply(CMLT[,9:13],1,sum)))
(PYt=cbind(apply(PYT[,1:4],1,sum),apply(PYT[,5:8],1,sum),apply(PYT[,9:13],1,sum)))
(AGEt=cbind(apply(AGET[,1:4],1,mean),apply(AGET[,5:8],1,mean),apply(AGET[,9:13],1,mean)))
(AGEw=cbind(apply(AGET[,1:4]*PYT[,1:4],1,sum),apply(AGET[,5:8]*PYT[,5:8],1,sum),apply(AGET[,9:13]*PYT[,9:13],1,sum)))
(AGEtw=AGEw/PYt)
(incid=CMLt/PYt)
#	(incidUpper=(CMLt+2*sqrt(CMLt))/PYt)
#	(incidLow=(CMLt-2*sqrt(CMLt))/PYt)
brb=c("blue","red","black")
pchs=c(2,1,22)
graphics.off()
if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
if(length(grep("darwin",R.Version()$os))) windows <- function( ... ) quartz( ... )

windows(width=6,height=6)
par(mfrow=c(1,1),mar=c(4.7,5.4,1.3,.8))
matplot(t(AGEtw),t(incid),log="y",type='b',ylab=expression(paste("Cases per ",10^5," Person-Years")),
        xlab="Age at exposure",
        ylim=c(1e-6,8e-4),lty=1,col=brb,cex=2.5,pch=pchs,lwd=3,cex.lab=2,cex.axis=2,yaxt="n")
axis(2,at=c(1e-6,1e-5,1e-4),labels=c(0.1,1,10),cex.axis=2)
text(30,18e-5,"mostly radiogenic",cex=1.5)
text(37,1.3e-5,"mostly background",cex=1.5,col="blue")
legend(20,0.6e-5,c("High Dose","Medium Dose","Low Dose"),pch=pchs[3:1],col=brb[3:1],cex=1.5,pt.cex=2,pt.lwd=2)

