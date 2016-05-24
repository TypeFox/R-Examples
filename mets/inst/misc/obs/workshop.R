
###############################
## installation, R (>=2.12.0)
###############################

install.packages("mets")

###############################
## Load simulated data
###############################

library(mets)
data(prt)

###############################
## Estimation of cumulative incidence
###############################

times <- seq(40,100,by=2)
cifmod <- comp.risk(Surv(time,status>0)~+1+cluster(id),data=prt,prt$status,causeS=2,n.sim=0,
                  times=times,conservative=1,max.clust=NULL,model="fg")

theta.des <- model.matrix(~-1+factor(zyg),data=prt) ## design for MZ/DZ status
or1 <- or.cif(cifmod,data=prt,cause1=2,cause2=2,theta.des=theta.des,
              score.method="fisher.scoring")
summary(or1)
or1$score

pcif <- predict(cifmod,X=1,resample.iid=0,uniform=0,se=0)

png(filename="pcif.png")
plot(pcif,multiple=1,se=0,uniform=0,ylim=c(0,0.15))
dev.off()

###############################
## Correcting for country
###############################

png(filename="pcifl.png")
table(prt$country)

times <- seq(40,100,by=2)
cifmodl <-comp.risk(Surv(time,status>0)~-1+factor(country)+cluster(id),data=prt,
                    prt$status,causeS=2,n.sim=0,times=times,conservative=1,
                    max.clust=NULL,cens.model="aalen")
pcifl <- predict(cifmodl,X=diag(4),se=0,uniform=0)
plot(pcifl,multiple=1,se=0,uniform=0,col=1:4,ylim=c(0,0.2))
legend("topleft",levels(prt$country),col=1:4,lty=1)
dev.off()

theta.des <- model.matrix(~-1+factor(zyg),data=prt) ## design for MZ/DZ status
or.country <- or.cif(cifmodl,data=prt,cause1=2,cause2=2,theta.des=theta.des,
                     theta=c(2.8,6.9),score.method="fisher.scoring")
summary(or.country)

cifmodlr <-comp.risk(Surv(time,status>0)~+1+const(factor(country))+cluster(id),data=prt,
                    prt$status,causeS=2,n.sim=0,times=times,conservative=1,max.clust=NULL,model="fg",
                    cens.model="aalen",cens.formula=~~factor(country))
pciflr <- predict(cifmodlr,X=rep(1,4),Z=rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(0,0,1)),se=0,uniform=0)

png(filename="pcif2.png")
par(mfrow=c(1,2))
plot(pcifl,multiple=1,se=0,uniform=0,col=1:4,ylim=c(0,0.2))
legend("topleft",levels(prt$country),col=1:4,lty=1)
plot(pciflr,multiple=1,se=0,uniform=0,col=1:4,ylim=c(0,0.2))
legend("topleft",levels(prt$country),col=1:4,lty=1)
dev.off()

or.countryr <- or.cif(cifmodlr,data=prt,cause1=2,cause2=2,theta.des=theta.des,
                     theta=c(2.8,6.9),score.method="fisher.scoring")
summary(or.countryr)

###############################
## Concordance estimation
###############################

### ignoring country 
p33 <- bicomprisk(Hist(time,status)~strata(zyg)+id(id),data=prt,cause=c(2,2),return.data=1)

p33dz <- p33$model$"DZ"$comp.risk
p33mz <- p33$model$"MZ"$comp.risk

png(filename="p33dz.png")
plot(p33dz,ylim=c(0,0.1))
dev.off()

png(filename="pcaconc.png")
par(new=TRUE)
plot(p33mz,ylim=c(0,0.1),col=3)
title(main="Concordance Prostate cancer")
lines(pcif$time,pcif$P1^2,col=2)
dev.off()

### test for genetic effect 
conc1 <- p33dz
conc2 <- p33mz

test.conc(p33dz,p33mz);

data33mz <- p33$model$"MZ"$data
data33mz$zyg <- 1
data33dz <- p33$model$"DZ"$data
data33dz$zyg <- 0
data33 <- rbind(data33mz,data33dz)

library(cmprsk)
ftime <- data33$time
fstatus <- data33$status
table(fstatus)

group <- data33$zyg
graytest <- cuminc(ftime,fstatus,group)
graytest

zygeffect <- comp.risk(Surv(time,status==0)~const(zyg),
                  data=data33,data33$status,causeS=1,
                  cens.model="aalen",model="logistic",conservative=1)
summary(zygeffect)

png(filename="casewise.png")
case33mz <- conc2case(p33mz,pcif)
case33dz <- conc2case(p33dz,pcif)

plot(case33mz$casewise,se=0,col=2)
par(new=TRUE)
plot(case33dz$casewise,se=0)
dev.off()

###############################
## Effect of zygosity correcting for country
###############################

p33l <- bicomprisk(Hist(time,status)~country+strata(zyg)+id(id),
                data=prt,cause=c(2,2),return.data=1,robust=1)

data33mz <- p33l$model$"MZ"$data
data33mz$zyg <- 1
data33dz <- p33l$model$"DZ"$data
data33dz$zyg <- 0
data33 <- rbind(data33mz,data33dz)

zygeffectl <- comp.risk(Surv(time,status==0)~const(country)+const(zyg),
                  data=data33,data33$status,causeS=1,
                  cens.model="aalen",model="logistic",conservative=1)
summary(zygeffectl)

zygeffectpl <- comp.risk(Surv(time,status==0)~const(country)+const(zyg),
                  data=data33,data33$status,causeS=1,
                  cens.model="aalen",model="fg",conservative=1)
summary(zygeffectpl)

zygeffectll <- comp.risk(Surv(time,status==0)~country+const(zyg),
                         data=data33,data33$status,causeS=1,
                         cens.model="aalen",model="logistic",conservative=1)
summary(zygeffectll)

###############################
## Liability model, ignoring censoring
###############################

(M <- with(prt, table(cancer,zyg)))

coef(lm(cancer~-1+zyg,prt))

## Saturated model
bpmz <- 
    biprobit(cancer~1 + cluster(id), 
             data=subset(prt,zyg=="MZ"), eqmarg=TRUE)

logLik(bpmz) # Log-likelihood
AIC(bpmz) # AIC
coef(bpmz) # Parameter estimates
vcov(bpmz) # Asymptotic covariance
summary(bpmz) # concordance, case-wise, tetrachoric correlations, ...

bp0 <- biprobit(cancer~1 + cluster(id)+strata(zyg), data=prt)

summary(bp0)

## Eq. marginals MZ/DZ
bp1 <- bptwin(cancer~1,zyg="zyg",DZ="DZ",id="id",type="u",data=prt)
summary(bp1) # Components (concordance,cor,...) can be extracted from returned list

compare(bp0,bp1) # LRT

## Polygenic model
args(bptwin)

bp2 <- bptwin(cancer~1,zyg="zyg",DZ="DZ",id="id",type="ace",data=prt)
summary(bp2)

###############################
## Liability model, IPCW
###############################

png(filename="ipw.png")
## Probability weights based on Aalen's additive model 
prtw <- ipw(Surv(time,status==0)~country, data=prt,
            cluster="id",weightname="w") 
plot(0,type="n",xlim=range(prtw$time),ylim=c(0,1),xlab="Age",ylab="Probability")
count <- 0
for (l in unique(prtw$country)) {
    count <- count+1
    prtw <- prtw[order(prtw$time),]
    with(subset(prtw,country==l), 
         lines(time,w,col=count,lwd=2))
}
legend("topright",legend=unique(prtw$country),col=1:4,pch=1)
dev.off()

bpmzIPW <- 
              biprobit(cancer~1 + cluster(id), 
                       data=subset(prtw,zyg=="MZ"), 
                       weight="w")
(smz <- summary(bpmzIPW))

png(filename="cif2.png")
## CIF
plot(pcif,multiple=1,se=0,uniform=0,ylim=c(0,0.15))
abline(h=smz$prob["Marginal",],lwd=c(2,1,1))
## Wrong estimates:
abline(h=summary(bpmz)$prob["Marginal",],lwd=c(2,1,1),col="lightgray")
dev.off()

png(filename="conc2.png")
## Concordance
plot(p33mz,ylim=c(0,0.1))
abline(h=smz$prob["Concordance",],lwd=c(2,1,1))
## Wrong estimates:
abline(h=summary(bpmz)$prob["Concordance",],lwd=c(2,1,1),col="lightgray")
dev.off()

bp3 <- bptwin(cancer~1,zyg="zyg",DZ="DZ",id="id",
              type="ace",data=prtw,weight="w")
summary(bp3)

bp4 <- bptwin(cancer~1,zyg="zyg",DZ="DZ",id="id",
              type="u",data=prtw,weight="w")
summary(bp4)

score(bp4) ## Check convergence

bp5 <- bptwin(cancer~1,zyg="zyg",DZ="DZ",id="id",
              type="ade",data=prtw,weight="w")
summary(bp5)

###############################
## Adjusting for covariates
###############################

bp6 <- bptwin(cancer~country,zyg="zyg",DZ="DZ",id="id",
              type="ace",data=prtw,weight="w")
summary(bp6)

bp7 <- bptwin(cancer~country,zyg="zyg",DZ="DZ",id="id",
              type="u",data=prtw,weight="w")
summary(bp7)

bp8 <- bptwin(cancer~strata(country),zyg="zyg",DZ="DZ",id="id",
              type="u",data=prtw,weight="w")

summary(bp8)

## Wald test
B <- (lava::contrmat(3,4))[-(1:3),]
compare(bp8,contrast=B)

###############################
## Cumulative heritability
###############################

args(cumh)

ch1 <- cumh(cancer~1,time="time",zyg="zyg",DZ="DZ",id="id",
            type="ace",data=prtw,weight="w")
#+BEGIN_SRC R
summary(ch1)

png(filename="cumh.png")
plot(ch1)
dev.off()
