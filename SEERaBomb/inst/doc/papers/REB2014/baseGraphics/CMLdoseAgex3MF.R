# Figure 3 of current MS. Differs from Fig 6 in Blood in that gender is separated and city can be too
rm(list=ls(all=TRUE))
cols=c("city","sex","doseg","agexg","calg","kerma","PY","adjPY","num.entering",
"age","agex","tsx","cal","sv","gam","neut","lymphoma","NHL","leukemia","AML","ALL","CML","ATL","MML")      
d<-read.table("/data/abomb/HEMA87.DAT", header=F,col.names=cols);
d=d[d$adjPY>0,] #remove any recs with zero py
d$adjPY=d$adjPY*1e4 # convert to real person years
d=d[d$kerma==1,] # take only kerma < 4 Gy
CITY=1  # if 0 pool both cites, else 1 = H only and 2 = N only
if (CITY>0) d=d[d$city==CITY,]#hiroshima = 1
m=d[d$sex==1,]; f=d[d$sex==2,] 
getMats<-function(m) {
	svc<-with(m,cut(sv,c(-1,.02,1,10)))
	PYT<-with(m,tapply(adjPY,list(svc,agexg),sum,na.rm=TRUE))
	AGET<-with(m,tapply(agex,list(svc,agexg),mean,na.rm=TRUE))
	CMLT<-with(m,tapply(CML,list(svc,agexg),sum,na.rm=TRUE)) 
	(CMLt=cbind(apply(CMLT[,1:4],1,sum),apply(CMLT[,5:8],1,sum),apply(CMLT[,9:13],1,sum)))
	(PYt=cbind(apply(PYT[,1:4],1,sum),apply(PYT[,5:8],1,sum),apply(PYT[,9:13],1,sum)))
	(AGEt=cbind(apply(AGET[,1:4],1,mean),apply(AGET[,5:8],1,mean),apply(AGET[,9:13],1,mean)))
	(AGEw=cbind(apply(AGET[,1:4]*PYT[,1:4],1,sum),apply(AGET[,5:8]*PYT[,5:8],1,sum),apply(AGET[,9:13]*PYT[,9:13],1,sum)))
	list(AGEtw=AGEw/PYt,incid=CMLt/PYt)}
M=getMats(m)
F=getMats(f)
brb=c("blue","red","black");pchs=c(2,1,22)
graphics.off()
if(length(grep("linux",R.Version()$os))) windows <- function( ... ) X11( ... )
windows(width=10,height=5)
par(mfrow=c(1,2),mar=c(4.7,0,0.5,0.9),lwd=3,cex.lab=1.8,cex.axis=2,cex.main=1.7,oma=c(0,5.5,1.4,0.7))
matplot(t(M$AGEtw),t(M$incid),log="y",type='b',ylab="",	xlab="age at exposure",
    ylim=c(1e-6,15e-4),lty=1,col=brb,cex=2.5,pch=pchs,lwd=3,cex.lab=2,yaxt="n")
axis(2,at=c(1e-6,1e-5,1e-4,1e-3),labels=c(0.1,1,10,100),las=1,cex.axis=1.7,outer=T)
mtext(expression(paste("Cases per ",10^5," Person-Year-Sv")),side=2,line=3.3,cex=1.8,outer=T)
mtext("Males",side=3,line=-2,cex=1.8,adj=0.05); text(30,25e-5,"mostly radiogenic",cex=1.5)
legend(20,1.3e-5,c("High Dose","Medium Dose","Low Dose"),pch=pchs[3:1],col=brb[3:1],cex=1.5,pt.cex=2,pt.lwd=2,bty="n")
title("Japanese A-bomb Survivors",outer=T)
par(mar=c(4.7,1,0.5,0))
matplot(t(F$AGEtw),t(F$incid),log="y",type='b',ylab="",	xlab="age at exposure",
		ylim=c(1e-6,15e-4),lty=1,col=brb,cex=2.5,pch=pchs,lwd=3,cex.lab=2,yaxt="n")
text(30,25e-5,"mostly radiogenic",cex=1.5); mtext("Females",side=3,line=-2,cex=1.8,adj=0.05)