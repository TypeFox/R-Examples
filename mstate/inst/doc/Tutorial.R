### R code from vignette source 'Tutorial.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: mstate
###################################################
library(mstate)


###################################################
### code chunk number 2: mstate
###################################################
R.version$version.string
packageDescription("mstate", fields = "Version")


###################################################
### code chunk number 3: data
###################################################
data(ebmt3)
head(ebmt3)


###################################################
### code chunk number 4: tables
###################################################
n <- nrow(ebmt3)
table(ebmt3$dissub); round(100*table(ebmt3$dissub)/n)


###################################################
### code chunk number 5: tables2 (eval = FALSE)
###################################################
## table(ebmt3$age); round(100*table(ebmt3$age)/n)
## table(ebmt3$drmatch); round(100*table(ebmt3$drmatch)/n)
## table(ebmt3$tcd); round(100*table(ebmt3$tcd)/n)


###################################################
### code chunk number 6: tmat
###################################################
tmat <- matrix(NA,3,3)
tmat[1,2:3] <- 1:2
tmat[2,3] <- 3
dimnames(tmat) <- list(from=c("Tx","PR","RelDeath"),to=c("Tx","PR","RelDeath"))
tmat


###################################################
### code chunk number 7: transMat
###################################################
tmat <- transMat(x=list(c(2,3),c(3),c()), names=c("Tx","PR","RelDeath"))
tmat


###################################################
### code chunk number 8: tmat2
###################################################
tmat <- trans.illdeath(names=c("Tx","PR","RelDeath"))
tmat


###################################################
### code chunk number 9: paths
###################################################
paths(tmat)


###################################################
### code chunk number 10: times
###################################################
ebmt3$prtime <- ebmt3$prtime/365.25
ebmt3$rfstime <- ebmt3$rfstime/365.25


###################################################
### code chunk number 11: msprep
###################################################
covs <- c("dissub", "age", "drmatch", "tcd", "prtime")
msbmt <- msprep(time=c(NA,"prtime","rfstime"),
    status=c(NA,"prstat","rfsstat"), data=ebmt3, trans=tmat,
    keep=covs)


###################################################
### code chunk number 12: msbmt
###################################################
head(msbmt)


###################################################
### code chunk number 13: events
###################################################
events(msbmt)


###################################################
### code chunk number 14: expandcovs
###################################################
expcovs <- expand.covs(msbmt,covs[2:3],append=FALSE)
head(expcovs)


###################################################
### code chunk number 15: expandcovs2
###################################################
msbmt <- expand.covs(msbmt,covs,append=TRUE,longnames=FALSE)
head(msbmt)


###################################################
### code chunk number 16: Markovnonprop
###################################################
c1 <- coxph(Surv(Tstart,Tstop,status) ~ dissub1.1 + dissub2.1 + age1.1 + age2.1 + drmatch.1 + tcd.1 
    + dissub1.2 + dissub2.2 + age1.2 + age2.2 + drmatch.2 + tcd.2 
    + dissub1.3 + dissub2.3 + age1.3 + age2.3 + drmatch.3 + tcd.3 
    + strata(trans),
    data=msbmt, method="breslow")
c1


###################################################
### code chunk number 17: Markovprop
###################################################
msbmt$pr <- 0
msbmt$pr[msbmt$trans==3] <- 1
c2 <- coxph(Surv(Tstart,Tstop,status) ~ dissub1.1 + dissub2.1 + age1.1 + age2.1 + drmatch.1 + tcd.1
    + dissub1.2 + dissub2.2 + age1.2 + age2.2 + drmatch.2 + tcd.2
    + dissub1.3 + dissub2.3 + age1.3 + age2.3 + drmatch.3 + tcd.3
    + pr
    + strata(to),
    data=msbmt, method="breslow")
c2


###################################################
### code chunk number 18: coxzph
###################################################
cox.zph(c2)


###################################################
### code chunk number 19: nonMarkovPH
###################################################
c3 <- coxph(Surv(Tstart,Tstop,status) ~ dissub1.1 + dissub2.1 + age1.1 + age2.1 + drmatch.1 + tcd.1
    + dissub1.2 + dissub2.2 + age1.2 + age2.2 + drmatch.2 + tcd.2
    + dissub1.3 + dissub2.3 + age1.3 + age2.3 + drmatch.3 + tcd.3
    + pr
    + prtime.3
    + strata(to),
    data=msbmt, method="breslow")
c3


###################################################
### code chunk number 20: semiMarkov
###################################################
c4 <- coxph(Surv(time,status) ~ dissub1.1 + dissub2.1 + age1.1 + age2.1 + drmatch.1 + tcd.1 
    + dissub1.2 + dissub2.2 + age1.2 + age2.2 + drmatch.2 + tcd.2 
    + dissub1.3 + dissub2.3 + age1.3 + age2.3 + drmatch.3 + tcd.3 
    + strata(trans),
    data=msbmt, method="breslow")
c5 <- coxph(Surv(time,status) ~ dissub1.1 + dissub2.1 + age1.1 + age2.1 + drmatch.1 + tcd.1
    + dissub1.2 + dissub2.2 + age1.2 + age2.2 + drmatch.2 + tcd.2
    + dissub1.3 + dissub2.3 + age1.3 + age2.3 + drmatch.3 + tcd.3
    + pr
    + strata(to),
    data=msbmt, method="breslow")
c6 <- coxph(Surv(time,status) ~ dissub1.1 + dissub2.1 + age1.1 + age2.1 + drmatch.1 + tcd.1
    + dissub1.2 + dissub2.2 + age1.2 + age2.2 + drmatch.2 + tcd.2
    + dissub1.3 + dissub2.3 + age1.3 + age2.3 + drmatch.3 + tcd.3
    + pr
    + prtime.3
    + strata(to),
    data=msbmt, method="breslow")


###################################################
### code chunk number 21: newdata
###################################################
newd <- data.frame(dissub=rep(0,3),age=rep(0,3),drmatch=rep(0,3),tcd=rep(0,3),trans=1:3)
newd$dissub <- factor(newd$dissub,levels=0:2,labels=levels(ebmt3$dissub))
newd$age <- factor(newd$age,levels=0:2,labels=levels(ebmt3$age))
newd$drmatch <- factor(newd$drmatch,levels=0:1,labels=levels(ebmt3$drmatch))
newd$tcd <- factor(newd$tcd,levels=0:1,labels=levels(ebmt3$tcd))
attr(newd, "trans") <- tmat
class(newd) <- c("msdata","data.frame")
newd <- expand.covs(newd,covs[1:4],longnames=FALSE)
newd$strata=1:3
newd


###################################################
### code chunk number 22: msfit
###################################################
msf1 <- msfit(c1, newdata=newd, trans=tmat)


###################################################
### code chunk number 23: msfitsummary
###################################################
summary(msf1)


###################################################
### code chunk number 24: msf1cov
###################################################
vH1 <- msf1$varHaz
head(vH1[vH1$trans1==1 & vH1$trans2==1,])
tail(vH1[vH1$trans1==1 & vH1$trans2==1,])
tail(vH1[vH1$trans1==1 & vH1$trans2==2,])
tail(vH1[vH1$trans1==1 & vH1$trans2==3,])
tail(vH1[vH1$trans1==2 & vH1$trans2==3,])


###################################################
### code chunk number 25: msfit2
###################################################
newd$strata=c(1,2,2)
newd$pr <- c(0,0,1)
msf2 <- msfit(c2, newdata=newd, trans=tmat)
summary(msf2)
vH2 <- msf2$varHaz
tail(vH2[vH2$trans1==1 & vH2$trans2==2,])
tail(vH2[vH2$trans1==1 & vH2$trans2==3,])
tail(vH2[vH2$trans1==2 & vH2$trans2==3,])


###################################################
### code chunk number 26: fig14a (eval = FALSE)
###################################################
## par(mfrow=c(1,2))
## plot(msf1,cols=rep(1,3),lwd=2,lty=1:3,
##     xlab="Years since transplant",ylab="Stratified baseline hazards",legend.pos=c(2,0.9))
## plot(msf2,cols=rep(1,3),lwd=2,lty=1:3,
##     xlab="Years since transplant",ylab="Proportional baseline hazards",legend.pos=c(2,0.9))
## par(mfrow=c(1,1))


###################################################
### code chunk number 27: fig14b
###################################################
par(mfrow=c(1,2))
plot(msf1,cols=rep(1,3),lwd=2,lty=1:3,
    xlab="Years since transplant",ylab="Stratified baseline hazards",legend.pos=c(2,0.9))
plot(msf2,cols=rep(1,3),lwd=2,lty=1:3,
    xlab="Years since transplant",ylab="Proportional baseline hazards",legend.pos=c(2,0.9))
par(mfrow=c(1,1))


###################################################
### code chunk number 28: probtrans1
###################################################
pt <- probtrans(msf2, predt=0)


###################################################
### code chunk number 29: pt13
###################################################
head(pt[[3]])
tail(pt[[3]])


###################################################
### code chunk number 30: pt12
###################################################
summary(pt, from=2)


###################################################
### code chunk number 31: pt11
###################################################
summary(pt, from=1)


###################################################
### code chunk number 32: tmat2
###################################################
tmat2 <- transMat(x=list(c(2,4),c(3),c(),c()))
tmat2


###################################################
### code chunk number 33: probtrans2
###################################################
msf2$trans <- tmat2
pt <- probtrans(msf2, predt=0)
summary(pt, from=1)


###################################################
### code chunk number 34: pt2plot (eval = FALSE)
###################################################
## plot(pt,ord=c(2,3,4,1),lwd=2,
##     xlab="Years since transplant",ylab="Prediction probabilities",cex=0.75,
##     legend=c("Alive in remission, no PR","Alive in remission, PR",
##     "Relapse or death after PR","Relapse or death without PR"))


###################################################
### code chunk number 35: pt2plot
###################################################
plot(pt,ord=c(2,3,4,1),lwd=2,
    xlab="Years since transplant",ylab="Prediction probabilities",cex=0.75,
    legend=c("Alive in remission, no PR","Alive in remission, PR",
    "Relapse or death after PR","Relapse or death without PR"))


###################################################
### code chunk number 36: probtrans2
###################################################
pt <- probtrans(msf2, predt=0.5)
summary(pt, from=1)


###################################################
### code chunk number 37: pt2plot2state (eval = FALSE)
###################################################
## plot(pt,ord=c(2,3,4,1),lwd=2,
##     xlab="Years since transplant",ylab="Prediction probabilities",cex=0.75,
##     legend=c("Alive in remission, no PR","Alive in remission, PR",
##     "Relapse or death after PR","Relapse or death without PR"))


###################################################
### code chunk number 38: pt2plot2do
###################################################
plot(pt,ord=c(2,3,4,1),lwd=2,
    xlab="Years since transplant",ylab="Prediction probabilities",cex=0.75,
    legend=c("Alive in remission, no PR","Alive in remission, PR",
    "Relapse or death after PR","Relapse or death without PR"))


###################################################
### code chunk number 39: makefig17a
###################################################
msf2$trans <- tmat
msf.20 <- msf2 # copy msfit result for reference (young) patient
newd <- newd[,1:5] # use the basic covariates of the reference patient
newd2 <- newd
newd2$age <- 1
newd2$age <- factor(newd2$age,levels=0:2,labels=levels(ebmt3$age))
attr(newd2, "trans") <- tmat
class(newd2) <- c("msdata","data.frame")
newd2 <- expand.covs(newd2,covs[1:4],longnames=FALSE)
newd2$strata=c(1,2,2)
newd2$pr <- c(0,0,1)
msf.2040 <- msfit(c2, newdata=newd2, trans=tmat)
newd3 <- newd
newd3$age <- 2
newd3$age <- factor(newd3$age,levels=0:2,labels=levels(ebmt3$age))
attr(newd3, "trans") <- tmat
class(newd3) <- c("msdata","data.frame")
newd3 <- expand.covs(newd3,covs[1:4],longnames=FALSE)
newd3$strata=c(1,2,2)
newd3$pr <- c(0,0,1)
msf.40 <- msfit(c2, newdata=newd3, trans=tmat)
pt.20 <- probtrans(msf.20,predt=0) # original young (<= 20) patient
pt.201 <- pt.20[[1]]; pt.202 <- pt.20[[2]]
pt.2040 <- probtrans(msf.2040,predt=0) # patient 20-40
pt.20401 <- pt.2040[[1]]; pt.20402 <- pt.2040[[2]]
pt.40 <- probtrans(msf.40,predt=0) # patient > 40
pt.401 <- pt.40[[1]]; pt.402 <- pt.40[[2]]


###################################################
### code chunk number 40: rfs5
###################################################
pt.201[488:489,] # 5 years falls between 488th and 489th time point
pt.202[488:489,] # 5-years probabilities


###################################################
### code chunk number 41: makefig17b (eval = FALSE)
###################################################
## plot(pt.201$time,1-pt.201$pstate3,ylim=c(0.425,1),type="s",lwd=2,col="red",xlab="Years since transplant",ylab="Relapse-free survival")
## lines(pt.20401$time,1-pt.20401$pstate3,type="s",lwd=2,col="blue")
## lines(pt.401$time,1-pt.401$pstate3,type="s",lwd=2,col="green")
## lines(pt.202$time,1-pt.202$pstate3,type="s",lwd=2,col="red",lty=2)
## lines(pt.20402$time,1-pt.20402$pstate3,type="s",lwd=2,col="blue",lty=2)
## lines(pt.402$time,1-pt.402$pstate3,type="s",lwd=2,col="green",lty=2)
## legend(6,1,c("no PR","PR"),lwd=2,lty=1:2,xjust=1,bty="n")
## legend("topright",c("<=20","20-40",">40"),lwd=2,col=c("red","blue","green"),bty="n")


###################################################
### code chunk number 42: makefig17c
###################################################
plot(pt.201$time,1-pt.201$pstate3,ylim=c(0.4,1),type="s",lwd=2,col="red",xlab="Years since transplant",ylab="Relapse-free survival")
lines(pt.20401$time,1-pt.20401$pstate3,type="s",lwd=2,col="blue")
lines(pt.401$time,1-pt.401$pstate3,type="s",lwd=2,col="green")
lines(pt.202$time,1-pt.202$pstate3,type="s",lwd=2,col="red",lty=2)
lines(pt.20402$time,1-pt.20402$pstate3,type="s",lwd=2,col="blue",lty=2)
lines(pt.402$time,1-pt.402$pstate3,type="s",lwd=2,col="green",lty=2)
legend(6,1,c("no PR","PR"),lwd=2,lty=1:2,xjust=1,bty="n")
legend("topright",c("<=20","20-40",">40"),lwd=2,col=c("red","blue","green"),bty="n")


###################################################
### code chunk number 43: backward
###################################################
pt.20 <- probtrans(msf.20,direction="fixedhorizon",predt=5)
pt.201 <- pt.20[[1]]; pt.202 <- pt.20[[2]]
head(pt.201); head(pt.202)


###################################################
### code chunk number 44: backward2
###################################################
pt.2040 <- probtrans(msf.2040,direction="fixedhorizon",predt=5) # patient 20-40
pt.20401 <- pt.2040[[1]]; pt.20402 <- pt.2040[[2]]
pt.40 <- probtrans(msf.40,direction="fixedhorizon",predt=5) # patient > 40
pt.401 <- pt.40[[1]]; pt.402 <- pt.40[[2]]


###################################################
### code chunk number 45: makenewa (eval = FALSE)
###################################################
## plot(pt.201$time,1-pt.201$pstate3,ylim=c(0.425,1),type="s",
##     lwd=2,col="red",xlab="Years since transplant",
##     ylab="Relapse-free survival")
## lines(pt.20401$time,1-pt.20401$pstate3,type="s",lwd=2,col="blue")
## lines(pt.401$time,1-pt.401$pstate3,type="s",lwd=2,col="green")
## lines(pt.202$time,1-pt.202$pstate3,type="s",lwd=2,col="red",lty=2)
## lines(pt.20402$time,1-pt.20402$pstate3,type="s",lwd=2,col="blue",lty=2)
## lines(pt.402$time,1-pt.402$pstate3,type="s",lwd=2,col="green",lty=2)
## legend("topleft",c("<=20","20-40",">40"),lwd=2,col=c("red","blue","green"),bty="n")
## legend(1,1,c("no PR","PR"),lwd=2,lty=1:2,bty="n")
## title(main="Backward prediction")


###################################################
### code chunk number 46: makenewb
###################################################
plot(pt.201$time,1-pt.201$pstate3,ylim=c(0.425,1),type="s",
    lwd=2,col="red",xlab="Prediction time",
    ylab="Relapse-free survival")
lines(pt.20401$time,1-pt.20401$pstate3,type="s",lwd=2,col="blue")
lines(pt.401$time,1-pt.401$pstate3,type="s",lwd=2,col="green")
lines(pt.202$time,1-pt.202$pstate3,type="s",lwd=2,col="red",lty=2)
lines(pt.20402$time,1-pt.20402$pstate3,type="s",lwd=2,col="blue",lty=2)
lines(pt.402$time,1-pt.402$pstate3,type="s",lwd=2,col="green",lty=2)
legend("topleft",c("<=20","20-40",">40"),lwd=2,col=c("red","blue","green"),bty="n")
legend(1,1,c("no PR","PR"),lwd=2,lty=1:2,bty="n")
title(main="Backward prediction")


###################################################
### code chunk number 47: aidssi
###################################################
data(aidssi)
si <- aidssi # Just a shorter name
head(si)
table(si$status)


###################################################
### code chunk number 48: tmatcr
###################################################
tmat <- trans.comprisk(2,names=c("event-free","AIDS","SI"))
tmat


###################################################
### code chunk number 49: dataprep
###################################################
si$stat1 <- as.numeric(si$status==1)
si$stat2 <- as.numeric(si$status==2)
silong <- msprep(time=c(NA,"time","time"),status=c(NA,"stat1","stat2"),data=si,
                    keep="ccr5",trans=tmat)


###################################################
### code chunk number 50: check
###################################################
events(silong)


###################################################
### code chunk number 51: dummies
###################################################
silong <- expand.covs(silong,"ccr5")
silong[1:8,] # shows the first four patients in long format, as in tutorial


###################################################
### code chunk number 52: naive
###################################################
c1 <- coxph(Surv(time,status) ~ 1, data=silong, subset = (trans==1), method="breslow")
c2 <- coxph(Surv(time,status) ~ 1, data=silong, subset = (trans==2), method="breslow")
h1 <- survfit(c1)
h1 <- data.frame(time=h1$time,surv=h1$surv)
h2 <- survfit(c2)
h2 <- data.frame(time=h2$time,surv=h2$surv)


###################################################
### code chunk number 53: fig2 (eval = FALSE)
###################################################
## idx1 <- (h1$time<13) # this restricts the plot to the first 13 years
## plot(c(0,h1$time[idx1],13),c(1,h1$surv[idx1],min(h1$surv[idx1])),type="s",
##     xlim=c(0,13),ylim=c(0,1),xlab="Years from HIV infection",ylab="Probability",lwd=2)
## idx2 <- (h2$time<13)
## lines(c(0,h2$time[idx2],13),c(0,1-h2$surv[idx2],max(1-h2$surv[idx2])),type="s",lwd=2)
## text(8,0.71,adj=0,"AIDS")
## text(8,0.32,adj=0,"SI")


###################################################
### code chunk number 54: fig2do
###################################################
idx1 <- (h1$time<13) # this restricts the plot to the first 13 years
plot(c(0,h1$time[idx1],13),c(1,h1$surv[idx1],min(h1$surv[idx1])),type="s",
    xlim=c(0,13),ylim=c(0,1),xlab="Years from HIV infection",ylab="Probability",lwd=2)
idx2 <- (h2$time<13)
lines(c(0,h2$time[idx2],13),c(0,1-h2$surv[idx2],max(1-h2$surv[idx2])),type="s",lwd=2)
text(8,0.71,adj=0,"AIDS")
text(8,0.32,adj=0,"SI")


###################################################
### code chunk number 55: Cuminc1
###################################################
ci <- Cuminc(time=si$time, status=si$status)


###################################################
### code chunk number 56: Cuminc2
###################################################
ci <- Cuminc(time="time", status="status", data=aidssi)


###################################################
### code chunk number 57: ci
###################################################
head(ci); tail(ci)


###################################################
### code chunk number 58: fig3 (eval = FALSE)
###################################################
## idx0 <- (ci$time<13)
## plot(c(0,ci$time[idx0],13),c(1,1-ci$CI.1[idx0],min(1-ci$CI.1[idx0])),type="s",
##   xlim=c(0,13),ylim=c(0,1),xlab="Years from HIV infection",ylab="Probability",lwd=2)
## idx1 <- (h1$time<13)
## lines(c(0,h1$time[idx1],13),c(1,h1$surv[idx1],min(h1$surv[idx1])),type="s",lwd=2,col=8)
## lines(c(0,ci$time[idx0],13),c(0,ci$CI.2[idx0],max(ci$CI.2[idx0])),type="s",lwd=2)
## idx2 <- (h2$time<13)
## lines(c(0,h2$time[idx2],13),c(0,1-h2$surv[idx2],max(1-h2$surv[idx2])),type="s",lwd=2,col=8)
## text(8,0.77,adj=0,"AIDS")
## text(8,0.275,adj=0,"SI")


###################################################
### code chunk number 59: fig3do
###################################################
idx0 <- (ci$time<13)
plot(c(0,ci$time[idx0],13),c(1,1-ci$CI.1[idx0],min(1-ci$CI.1[idx0])),type="s",
  xlim=c(0,13),ylim=c(0,1),xlab="Years from HIV infection",ylab="Probability",lwd=2)
idx1 <- (h1$time<13)
lines(c(0,h1$time[idx1],13),c(1,h1$surv[idx1],min(h1$surv[idx1])),type="s",lwd=2,col=8)
lines(c(0,ci$time[idx0],13),c(0,ci$CI.2[idx0],max(ci$CI.2[idx0])),type="s",lwd=2)
idx2 <- (h2$time<13)
lines(c(0,h2$time[idx2],13),c(0,1-h2$surv[idx2],max(1-h2$surv[idx2])),type="s",lwd=2,col=8)
text(8,0.77,adj=0,"AIDS")
text(8,0.275,adj=0,"SI")


###################################################
### code chunk number 60: fig4 (eval = FALSE)
###################################################
## idx0 <- (ci$time<13)
## plot(c(0,ci$time[idx0]),c(0,ci$CI.1[idx0]),type="s",
##   xlim=c(0,13),ylim=c(0,1),xlab="Years from HIV infection",ylab="Probability",lwd=2)
## lines(c(0,ci$time[idx0]),c(0,ci$CI.1[idx0]+ci$CI.2[idx0]),type="s",lwd=2)
## text(13,0.5*max(ci$CI.1[idx0]),adj=1,"AIDS")
## text(13,max(ci$CI.1[idx0])+0.5*max(ci$CI.2[idx0]),adj=1,"SI")
## text(13,0.5+0.5*max(ci$CI.1[idx0])+0.5*max(ci$CI.2[idx0]),adj=1,"Event-free")


###################################################
### code chunk number 61: fig4do
###################################################
idx0 <- (ci$time<13)
plot(c(0,ci$time[idx0]),c(0,ci$CI.1[idx0]),type="s",
  xlim=c(0,13),ylim=c(0,1),xlab="Years from HIV infection",ylab="Probability",lwd=2)
lines(c(0,ci$time[idx0]),c(0,ci$CI.1[idx0]+ci$CI.2[idx0]),type="s",lwd=2)
text(13,0.5*max(ci$CI.1[idx0]),adj=1,"AIDS")
text(13,max(ci$CI.1[idx0])+0.5*max(ci$CI.2[idx0]),adj=1,"SI")
text(13,0.5+0.5*max(ci$CI.1[idx0])+0.5*max(ci$CI.2[idx0]),adj=1,"Event-free")


###################################################
### code chunk number 62: csh1
###################################################
coxph(Surv(time, status == 1) ~ ccr5, data = si) # AIDS
coxph(Surv(time, status == 2) ~ ccr5, data = si) # SI appearance


###################################################
### code chunk number 63: csh2
###################################################
coxph(Surv(time,status) ~ ccr5, data=silong, subset = (trans==1), method="breslow")
coxph(Surv(time,status) ~ ccr5, data=silong, subset = (trans==2), method="breslow")


###################################################
### code chunk number 64: csh3
###################################################
coxph(Surv(time, status) ~ ccr5WM.1 + ccr5WM.2 + strata(trans), data = silong)


###################################################
### code chunk number 65: csh4
###################################################
coxph(Surv(time, status) ~ ccr5 * factor(trans) + strata(trans), data = silong)


###################################################
### code chunk number 66: csh5
###################################################
coxph(Surv(time, status) ~ ccr5 + strata(trans), data = silong)


###################################################
### code chunk number 67: csh6
###################################################
coxph(Surv(time, status) ~ ccr5, data = silong)


###################################################
### code chunk number 68: csh7
###################################################
coxph(Surv(time, status != 0) ~ ccr5, data = si)


###################################################
### code chunk number 69: csh7
###################################################
coxph(Surv(time, status) ~ ccr5WM.1 + ccr5WM.2 + factor(trans), data = silong)


###################################################
### code chunk number 70: csh8
###################################################
coxph(Surv(time, status) ~ ccr5 * factor(trans), data = silong)


###################################################
### code chunk number 71: csh9
###################################################
coxph(Surv(time, status) ~ ccr5 * factor(trans) + cluster(id), data = silong)


###################################################
### code chunk number 72: ci1
###################################################
c1 <- coxph(Surv(time, status) ~ ccr5WM.1 + ccr5WM.2 + strata(trans), data = silong, method="breslow")


###################################################
### code chunk number 73: ci12
###################################################
WW <- data.frame(ccr5WM.1=c(0,0),ccr5WM.2=c(0,0),trans=c(1,2),strata=c(1,2))
msf.WW <- msfit(c1, WW, trans=tmat)


###################################################
### code chunk number 74: ci13
###################################################
pt.WW <- probtrans(msf.WW, 0)[[1]]


###################################################
### code chunk number 75: ci22
###################################################
WM <- data.frame(ccr5WM.1=c(1,0),ccr5WM.2=c(0,1),trans=c(1,2),strata=c(1,2))
msf.WM <- msfit(c1, WM, trans=tmat)
pt.WM <- probtrans(msf.WM, 0)[[1]]


###################################################
### code chunk number 76: fig5 (eval = FALSE)
###################################################
## idx1 <- (pt.WW$time<13)
## idx2 <- (pt.WM$time<13)
## ## AIDS
## plot(c(0,pt.WW$time[idx1]),c(0,pt.WW$pstate2[idx1]),type="s",ylim=c(0,0.5),xlab="Years from HIV infection",ylab="Probability",lwd=2)
## lines(c(0,pt.WM$time[idx2]),c(0,pt.WM$pstate2[idx2]),type="s",lwd=2,col=8)
## title(main="AIDS")
## text(9.2,0.345,"WW",adj=0,cex=0.75)
## text(9.2,0.125,"WM",adj=0,cex=0.75)
## ## SI appearance
## plot(c(0,pt.WW$time[idx1]),c(0,pt.WW$pstate3[idx1]),type="s",ylim=c(0,0.5),xlab="Years from HIV infection",ylab="Probability",lwd=2)
## lines(c(0,pt.WM$time[idx2]),c(0,pt.WM$pstate3[idx2]),type="s",lwd=2,col=8)
## title(main="SI appearance")
## text(7.5,0.31,"WW",adj=0,cex=0.75)
## text(7.5,0.245,"WM",adj=0,cex=0.75)


###################################################
### code chunk number 77: fig5do
###################################################
idx1 <- (pt.WW$time<13)
idx2 <- (pt.WM$time<13)
## AIDS
plot(c(0,pt.WW$time[idx1]),c(0,pt.WW$pstate2[idx1]),type="s",ylim=c(0,0.5),xlab="Years from HIV infection",ylab="Probability",lwd=2)
lines(c(0,pt.WM$time[idx2]),c(0,pt.WM$pstate2[idx2]),type="s",lwd=2,col=8)
title(main="AIDS")
text(9.2,0.345,"WW",adj=0,cex=0.75)
text(9.2,0.125,"WM",adj=0,cex=0.75)


###################################################
### code chunk number 78: fig5do2
###################################################
## SI appearance
plot(c(0,pt.WW$time[idx1]),c(0,pt.WW$pstate3[idx1]),type="s",ylim=c(0,0.5),xlab="Years from HIV infection",ylab="Probability",lwd=2)
lines(c(0,pt.WM$time[idx2]),c(0,pt.WM$pstate3[idx2]),type="s",lwd=2,col=8)
title(main="SI appearance")
text(7.5,0.31,"WW",adj=0,cex=0.75)
text(7.5,0.245,"WM",adj=0,cex=0.75)


###################################################
### code chunk number 79: fig7 (eval = FALSE)
###################################################
## ffs <- c(0,0.5,1,1.5,2,4)
## newmsf.WW <- msf.WW
## newmsf.WM <- msf.WM
## par(mfrow=c(2,3))
## for (ff in ffs) {
##     # WW
##     newmsf.WW$Haz$Haz[newmsf.WW$Haz$trans==1] <- ff*msf.WW$Haz$Haz[msf.WW$Haz$trans==1]
##     pt.WW <- probtrans(newmsf.WW, 0, variance=FALSE)[[1]]
##     # WM
##     newmsf.WM$Haz$Haz[newmsf.WM$Haz$trans==1] <- ff*msf.WM$Haz$Haz[msf.WM$Haz$trans==1]
##     pt.WM <- probtrans(newmsf.WM, 0, variance=FALSE)[[1]]
##     # Plot
##     idx1 <- (pt.WW$time<13)
##     idx2 <- (pt.WM$time<13)
##     plot(c(0,pt.WW$time[idx1]),c(0,pt.WW$pstate3[idx1]),type="s",ylim=c(0,0.52),xlab="Years from HIV infection",ylab="Probability",lwd=2)
##     lines(c(0,pt.WM$time[idx2]),c(0,pt.WM$pstate3[idx2]),type="s",lwd=2,col=8)
##     title(main=paste("Factor =",ff))
## }
## par(mfrow=c(1,1))


###################################################
### code chunk number 80: fig7do
###################################################
### Multiply baseline hazard of AIDS with factors (ff)
### of ff = 0, 0.5, 1, 1.5, 2, 4
ffs <- c(0,0.5,1,1.5,2,4)
newmsf.WW <- msf.WW
newmsf.WM <- msf.WM
par(mfrow=c(2,3))
for (ff in ffs) {
    # WW
    newmsf.WW$Haz$Haz[newmsf.WW$Haz$trans==1] <- ff*msf.WW$Haz$Haz[msf.WW$Haz$trans==1]
    pt.WW <- probtrans(newmsf.WW, 0, variance=FALSE)[[1]]
    # WM
    newmsf.WM$Haz$Haz[newmsf.WM$Haz$trans==1] <- ff*msf.WM$Haz$Haz[msf.WM$Haz$trans==1]
    pt.WM <- probtrans(newmsf.WM, 0, variance=FALSE)[[1]]
    # Plot
    idx1 <- (pt.WW$time<13)
    idx2 <- (pt.WM$time<13)
    plot(c(0,pt.WW$time[idx1]),c(0,pt.WW$pstate3[idx1]),type="s",ylim=c(0,0.52),xlab="Years from HIV infection",ylab="Probability",lwd=2)
    lines(c(0,pt.WM$time[idx2]),c(0,pt.WM$pstate3[idx2]),type="s",lwd=2,col=8)
    title(main=paste("Factor =",ff))
}
par(mfrow=c(1,1))


###################################################
### code chunk number 81: FG
###################################################
library(cmprsk)
sic <- si[!is.na(si$ccr5),]
ftime <- sic$time
fstatus <- sic$status
cov <- as.numeric(sic$ccr5)-1
# for failures of type 1 (AIDS)
z1 <- crr(ftime,fstatus,cov)
z1
# for failures of type 2 (SI)
z2 <- crr(ftime,fstatus,cov,failcode=2)
z2


###################################################
### code chunk number 82: fig8 (eval = FALSE)
###################################################
## z1.pr <- predict(z1,matrix(c(0,1),2,1))
## # this will contain predicted cum inc curves, both for WW (2nd column) and WM (3rd)
## z2.pr <- predict(z2,matrix(c(0,1),2,1))
## # Standard plots, not shown
## par(mfrow=c(1,2))
## plot(z1.pr,lty=1,lwd=2,color=c(8,1))
## plot(z2.pr,lty=1,lwd=2,color=c(8,1))
## par(mfrow=c(1,1))
## ## AIDS
## n1 <- nrow(z1.pr) # remove last jump
## plot(c(0,z1.pr[-n1,1]),c(0,z1.pr[-n1,2]),type="s",ylim=c(0,0.5),
##     xlab="Years from HIV infection",ylab="Probability",lwd=2)
## lines(c(0,z1.pr[-n1,1]),c(0,z1.pr[-n1,3]),type="s",lwd=2,col=8)
## title(main="AIDS")
## text(9.3,0.35,"WW",adj=0,cex=0.75)
## text(9.3,0.14,"WM",adj=0,cex=0.75)
## ## SI appearance
## n2 <- nrow(z2.pr) # again remove last jump
## plot(c(0,z2.pr[-n2,1]),c(0,z2.pr[-n2,2]),type="s",ylim=c(0,0.5),
##     xlab="Years from HIV infection",ylab="Probability",lwd=2)
## lines(c(0,z2.pr[-n2,1]),c(0,z2.pr[-n2,3]),type="s",lwd=2,col=8)
## title(main="SI appearance")
## text(7.9,0.28,"WW",adj=0,cex=0.75)
## text(7.9,0.31,"WM",adj=0,cex=0.75)


###################################################
### code chunk number 83: fig8do
###################################################
z1.pr <- predict(z1,matrix(c(0,1),2,1)) # this will contain predicted cum inc curves, both for WW (2nd column) and WM (3rd)
z2.pr <- predict(z2,matrix(c(0,1),2,1))
## AIDS
n1 <- nrow(z1.pr) # remove last jump
plot(c(0,z1.pr[-n1,1]),c(0,z1.pr[-n1,2]),type="s",ylim=c(0,0.5),xlab="Years from HIV infection",ylab="Probability",lwd=2)
lines(c(0,z1.pr[-n1,1]),c(0,z1.pr[-n1,3]),type="s",lwd=2,col=8)
title(main="AIDS")
text(9.3,0.35,"WW",adj=0,cex=0.75)
text(9.3,0.14,"WM",adj=0,cex=0.75)


###################################################
### code chunk number 84: fig8do2
###################################################
## SI appearance
n2 <- nrow(z2.pr) # again remove last jump
plot(c(0,z2.pr[-n2,1]),c(0,z2.pr[-n2,2]),type="s",ylim=c(0,0.5),xlab="Years from HIV infection",ylab="Probability",lwd=2)
lines(c(0,z2.pr[-n2,1]),c(0,z2.pr[-n2,3]),type="s",lwd=2,col=8)
title(main="SI appearance")
text(7.9,0.28,"WW",adj=0,cex=0.75)
text(7.9,0.31,"WM",adj=0,cex=0.75)


###################################################
### code chunk number 85: nonp
###################################################
ci <- Cuminc(si$time,si$status,group=si$ccr5)
ci.WW <- ci[ci$group=="WW",]
ci.WM <- ci[ci$group=="WM",]


###################################################
### code chunk number 86: Fig9 (eval = FALSE)
###################################################
## idx1 <- (ci.WW$time<13)
## idx2 <- (ci.WM$time<13)
## # AIDS
## plot(c(0,ci.WW$time[idx1]),c(0,ci.WW$CI.1[idx1]),type="s",ylim=c(0,0.5),xlab="Years from HIV infection",ylab="Probability",lwd=2)
## lines(c(0,ci.WM$time[idx2]),c(0,ci.WM$CI.1[idx2]),type="s",lwd=2,col=8)
## title(main="AIDS")
## text(9.3,0.35,"WW",adj=0,cex=0.75)
## text(9.3,0.11,"WM",adj=0,cex=0.75)
## # SI appearance
## plot(c(0,ci.WW$time[idx1]),c(0,ci.WW$CI.2[idx1]),type="s",ylim=c(0,0.5),xlab="Years from HIV infection",ylab="Probability",lwd=2)
## lines(c(0,ci.WM$time[idx2]),c(0,ci.WM$CI.2[idx2]),type="s",lwd=2,col=8)
## title(main="SI appearance")
## text(7.9,0.32,"WW",adj=0,cex=0.75)
## text(7.9,0.245,"WM",adj=0,cex=0.75)


###################################################
### code chunk number 87: Fig9do
###################################################
idx1 <- (ci.WW$time<13)
idx2 <- (ci.WM$time<13)
# AIDS
plot(c(0,ci.WW$time[idx1]),c(0,ci.WW$CI.1[idx1]),type="s",ylim=c(0,0.5),xlab="Years from HIV infection",ylab="Probability",lwd=2)
lines(c(0,ci.WM$time[idx2]),c(0,ci.WM$CI.1[idx2]),type="s",lwd=2,col=8)
title(main="AIDS")
text(9.3,0.35,"WW",adj=0,cex=0.75)
text(9.3,0.11,"WM",adj=0,cex=0.75)


###################################################
### code chunk number 88: Fig9do2
###################################################
# SI appearance
plot(c(0,ci.WW$time[idx1]),c(0,ci.WW$CI.2[idx1]),type="s",ylim=c(0,0.5),xlab="Years from HIV infection",ylab="Probability",lwd=2)
lines(c(0,ci.WM$time[idx2]),c(0,ci.WM$CI.2[idx2]),type="s",lwd=2,col=8)
title(main="SI appearance")
text(7.9,0.32,"WW",adj=0,cex=0.75)
text(7.9,0.245,"WM",adj=0,cex=0.75)


