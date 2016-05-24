### R code from vignette source 'hbimdetails.Rnw'

###################################################
### code chunk number 1: PreliminaryCalculations
###################################################
library(mvtnorm)
library(hbim)
NSIM<-5*10^5
SIGMAS.POWER<-c(9,65,5000)
SIGMAS<-log10(SIGMAS.POWER)/(2*qnorm(.975))
SCOLORS<-c("green","blue","red")
FACTORS<-c(1/10, 1/3, 1/2, 1)
FCOLORS<-c("red", "green", "blue", "black")
RHOS<-c(-.5,-.25,0, 0.25, 0.5, 0.75, 1)
RCOLORS<- c("purple","blue","black","green","yellow","red","black")
RECALCULATE.OUTPUT<-FALSE
#RECALCULATE.OUTPUT<-TRUE
time0<-proc.time()
if (RECALCULATE.OUTPUT){
    set.seed(1234521)
    MU<-((-40:40)/10)
    deff.sigma<-eff.sigma(mu=MU, sigmas=SIGMAS, COLORS = SCOLORS, rho = 0)
    deff.mu<-eff.mu(mu=MU, factor = FACTORS, COLORS = FCOLORS, sigma = SIGMAS[2], rho = 0)
    deff.rho<-eff.rho(mu=MU, sigma = SIGMAS[2], rho = RHOS, COLORS =RCOLORS,simulate=TRUE,nsim=NSIM)
    set.seed(32401)
    dpp.sigma<-pp.sigma(MU,sigmas=SIGMAS,COLORS = SCOLORS, rho = 0,nsim=NSIM)
    set.seed(21345123)
    dpp.mu<-pp.mu(MU,factor = FACTORS, COLORS = FCOLORS, sigma = SIGMAS[2], rho = 0, nsim=NSIM)
    set.seed(435919)
    dpp.rho<-pp.rho(MU,sigma = SIGMAS[2], rho = RHOS, COLORS =RCOLORS,nsim=NSIM)
     
}else{
    ## load from previously calculated
    data(deff.sigma)
    data(deff.mu)
    data(deff.rho)
    data(dpp.sigma)
    data(dpp.mu)
    data(dpp.rho)
}
time1<-proc.time()


###################################################
### code chunk number 2: hbimdetails.Rnw:355-357
###################################################
data(irdata)
irdata[35,]


###################################################
### code chunk number 3: hbimdetails.Rnw:367-370
###################################################
refnum<- irdata[35,"Reference"]
data(refs)
refs[refnum]


###################################################
### code chunk number 4: hbimdetails.Rnw:383-403
###################################################
calc.foldrange<-function(n,lower,upper,conf.level=.95){
    alpha<-1-conf.level
    s.t<- 
(sqrt(n)*(log10(upper)-log10(lower)))/(2*qt(1-alpha/2,n-1))
    s.z<-(sqrt(n)*(log10(upper)-log10(lower)))/(2*qnorm(1-alpha/2))
    range.t<- 10^(( qnorm(1-alpha/2) - qnorm(alpha/2) )*s.t )
    range.z<- 10^(( qnorm(1-alpha/2) - qnorm(alpha/2) )*s.z )
    col.names<-c("n",paste("lower",100*conf.level,"% cl"),
           paste("upper",100*conf.level,"% cl"),
           "standard deviation (log10 scale) by t",
           "standard deviation (log10 scale) by Z",
           paste(100*conf.level,"% fold-range, s estimated by t"),
           paste(100*conf.level,"% fold-range, s estimated by Z"))
    
out<-data.frame(n=n,lower=lower,upper=upper,s.byt=s.t,s.byz=s.z,
        foldrange.byt=range.t,foldrange.byz=range.z)

    #out<-list(out,row.names)
    out
}


###################################################
### code chunk number 5: hbimdetails.Rnw:406-407
###################################################
calc.foldrange(203,65.2,87.6)


###################################################
### code chunk number 6: hbimdetails.Rnw:412-414
###################################################
frall<-calc.foldrange(irdata$n,irdata$GMT.95.pct.interval.low.limit,irdata$GMT.95.pct.interval.high.limit)
frall[1:5,]


###################################################
### code chunk number 7: hbimdetails.Rnw:421-423
###################################################
quantile(frall[,"foldrange.byt"],probs=c(0,.025,.50,.975,1))
quantile(frall[,"foldrange.byz"],probs=c(0,.025,.50,.975,1))


###################################################
### code chunk number 8: hbimdetails.Rnw:426-428
###################################################
quantile(frall[,"s.byt"],probs=c(0,.025,.50,.975,1))
quantile(frall[,"s.byz"],probs=c(0,.025,.50,.975,1))


###################################################
### code chunk number 9: fig1
###################################################
do.fig1<-function(IRDATA){
    age<-as.character(IRDATA$Age.in.yrs.at.first.vaccination)
    antigen<-as.character(IRDATA$Antigen)
    uage<-sort(unique(age))
    pick.fig1<- age<"1" & antigen!="A" & antigen!="C" & antigen!="FIM" &
        antigen!="W135" & antigen!="Y" & antigen!="PRP*"

    alevels<-c("1","3","4","5","6B","7F","9V","14","18C","19F","23F","MenC",
        "PRP","DT","FHA","PRN","PT","TT","HBs","Polio-1","Polio-2","Polio-3")
    fig1.data<-IRDATA[pick.fig1,c("Antigen","n",
        "GMT.95.pct.interval.low.limit","GMT.95.pct.interval.high.limit")]
    cf<-calc.foldrange(fig1.data$n,fig1.data$GMT.95.pct.interval.low.limit,
        fig1.data$GMT.95.pct.interval.high.limit)

    nfig1<-dim(fig1.data)[[1]]
    x<-COL<-rep(NA,nfig1)
    bit<-.4
    bit2<-1
    space<-1
    xlevels<-c((0:10)*bit,10*bit+space+bit2*(0:7),10*bit+space+bit2*7+space+bit*(0:2))
    collevels<-c(rep("blue",11),"slateblue4","aquamarine","green","red","red1","red2","orange","black","pink3","pink","pink1")
    for (i in 1:length(alevels)){
        x[fig1.data$Antigen==alevels[i]]<-xlevels[i]
        COL[fig1.data$Antigen==alevels[i]]<-collevels[i]
    }

    fig1.data<-data.frame(antigen=ordered(fig1.data$Antigen,levels=alevels),
        foldrange=cf[,"foldrange.byt"]) 
    par(oma=c(4,0,0,0))
    plot(range(x),log10(range(fig1.data$foldrange)),type="n",ylab="95% fold range",xlab="",axes=FALSE)
    box()
    axis(2,at=c(1,2,3,4),labels=c("10","100","1000","10000"))
    axis(1,at=xlevels[c(6,12:19,21)],
        labels=c("Pneumonococcal","Meningococcal","HiB",
        "Diptheria Toxin","Pertussis FHA","Pertussis PRN",
        "Pertussis Toxin","Tetanus Toxin","HBSAg","Polio"),las=2,cex=.5)
    points(x,log10(fig1.data$foldrange),col=COL,cex=1.5,pch=18)
}
do.fig1(irdata)


###################################################
### code chunk number 10: hbimdetails.Rnw:500-501
###################################################
sigmas<-SIGMAS


###################################################
### code chunk number 11: hbimdetails.Rnw:516-520
###################################################
V<-matrix( c(sigmas[2],0,0,sigmas[2]),2,2) 
V
mu<-c(.1,.1)
hbrr(mu,V)


###################################################
### code chunk number 12: hbimdetails.Rnw:530-531
###################################################
SIGMAS


###################################################
### code chunk number 13: hbimdetails.Rnw:535-537 (eval = FALSE)
###################################################
## data(deff.sigma)
## plotlogm.resp(deff.sigma)


###################################################
### code chunk number 14: fig2
###################################################
plotlogm.resp(deff.sigma)


###################################################
### code chunk number 15: fig3a
###################################################
plotresp.mix(deff.sigma)


###################################################
### code chunk number 16: equiv.increase
###################################################
D<-deff.sigma
equiv.increase(D$mu,D$out1[,2],D$mu,D$out2[,2],.5)
D<-dpp.sigma
equiv.increase(D$mu,D$out1[,2],D$mu,D$out2[,2],50)


###################################################
### code chunk number 17: explanationfig3b
###################################################
par(mfrow=c(1,1))
plot(c(-2,2),c(0,1),axes=FALSE,xlab="Standardized Geometric Mean Antibody",ylab="Efficacy",type="n")
#axis(1,at=c(-2,-1,0,1,2),labels=as.character(c(.01,.1,1,10,100)))
#axis(2)
box()
lines(deff.sigma$mu,deff.sigma$out1[,1],lty=1,col=deff.sigma$col1[1])
lines(deff.sigma$mu,deff.sigma$out2[,1],lty=5,col=deff.sigma$col2[1])

a1<-0
e1<-deff.sigma$out1[,1][deff.sigma$mu==a1]
e2<-deff.sigma$out2[,1][deff.sigma$mu==a1]
NOTNA<-!is.na( deff.sigma$out1[,1] ) & !is.na(deff.sigma$mu )
a2<-approx(deff.sigma$out1[NOTNA,1],deff.sigma$mu[NOTNA],xout=e2)$y

lines(c(a1,a1),c(-1,e2),lty=2,col="grey")
lines(c(-4,a2),c(e2,e2),lty=2,col="grey")
lines(c(-4,a1),c(e1,e1),lty=2,col="grey")
lines(c(a2,a2),c(-1,e2),lty=2,col="grey")

axis(2,at=c(e1),labels=expression(e[1]))
axis(2,at=c(e2),labels=expression(e[2]))
axis(1,at=c(a1),labels=expression(a[1]))
axis(1,at=c(a2),labels=expression(a[2]))


###################################################
### code chunk number 18: fig3b
###################################################
plotresp.equiv(deff.sigma)


###################################################
### code chunk number 19: figppsigma
###################################################
plotlogm.resp(dpp.sigma,YLIM=c(0,100),YLAB="Percent Protected")


###################################################
### code chunk number 20: fig3c
###################################################
plotresp.mix(dpp.sigma,XYLIM=c(0,100),RLAB="% Protected by")


###################################################
### code chunk number 21: figfirsttry3d
###################################################
plotresp.equiv(dpp.sigma,XLIM=c(0,100),RLAB="% Protected by")


###################################################
### code chunk number 22: fig3d
###################################################
plotresp.equiv(dpp.sigma,XLIM=c(0,100),RLAB="% Protected by",bounds=c(.01,99.9))


###################################################
### code chunk number 23: hbimdetails.Rnw:715-716
###################################################
FACTORS


###################################################
### code chunk number 24: figLogmeanEffMu
###################################################
plotlogm.resp(deff.mu)


###################################################
### code chunk number 25: figLogmeanPPMu
###################################################
plotlogm.resp(dpp.mu,YLIM=c(0,100),YLAB="Percent Protected")


###################################################
### code chunk number 26: fig4a
###################################################
plotresp.mix(deff.mu)


###################################################
### code chunk number 27: fig4b
###################################################
plotresp.equiv(deff.mu)


###################################################
### code chunk number 28: fig4c
###################################################
plotresp.mix(dpp.mu,XYLIM=c(0,100),RLAB="% Protected by")


###################################################
### code chunk number 29: fig4d
###################################################
plotresp.equiv(dpp.mu,XLIM=c(0,100),RLAB="% Protected by",bounds=c(.01,99.9))


###################################################
### code chunk number 30: hbimdetails.Rnw:794-795
###################################################
RHOS


###################################################
### code chunk number 31: figLogmeanEffMu
###################################################
plotlogm.resp(deff.rho)


###################################################
### code chunk number 32: figLogmeanPPMu
###################################################
plotlogm.resp(dpp.rho,YLIM=c(0,100),YLAB="Percent Protected")


###################################################
### code chunk number 33: fig5a
###################################################
plotresp.mix(deff.rho)


###################################################
### code chunk number 34: fig5b
###################################################
plotresp.equiv(deff.rho)


###################################################
### code chunk number 35: fig5c
###################################################
plotresp.mix(dpp.rho,XYLIM=c(0,100),RLAB="% Protected by")


###################################################
### code chunk number 36: fig5d
###################################################
plotresp.equiv(dpp.rho,XLIM=c(0,100),RLAB="% Protected by",bounds=c(.01,99.9))


