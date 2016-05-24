simPV <-
function(se, sp, pr, n1, n0, NPV0, PPV0, conf.level=0.95, NSIM=500)
{

PPVt<-ppv(p=pr, se=se, sp=sp)
NPVt<-npv(p=pr, se=se, sp=sp)

x11<-rbinom(n=NSIM, size=n1, prob=se)
x01<-n1-x11
x00<-rbinom(n=NSIM, size=n0, prob=sp)
x10<-n0-x00

# x1 = c(x11, x01)
# x0 = c(x10, x00)

RDAT<-cbind(x11, x01, x10, x00)
I1<-1:2
I0<-3:4

RCINPV<-apply(X=RDAT, MARGIN=1,
 FUN=function(x){CombCInpv(x0=x[I0], x1=x[I1], p=pr, conf.level=conf.level, alternative="greater")$conf.int[1]})

NPVpower<-sum(RCINPV>NPV0)/NSIM
NPVcover<-sum(RCINPV<=NPVt)/NSIM
NPVlimit<-quantile(RCINPV, probs=c(0.1,0.2,0.5))

RCIPPV<-apply(X=RDAT, MARGIN=1,
 FUN=function(x){CombCIppv(x0=x[I0], x1=x[I1], p=pr, conf.level=conf.level, alternative="greater")$conf.int[1]})

PPVpower<-sum(RCIPPV>PPV0)/NSIM
PPVcover<-sum(RCIPPV<=PPVt)/NSIM
PPVlimit<-quantile(RCIPPV, probs=c(0.1,0.2,0.5))

OUTNPV<-c(NPVpower,NPVcover,NPVlimit, NPVt, NPV0)
OUTPPV<-c(PPVpower,PPVcover,PPVlimit, PPVt, PPV0)

OUT<-rbind(OUTNPV, OUTPPV)

rownames(OUT)<-c("NPV","PPV")
colnames(OUT)<-c("pow","cov","q10","q20","q50","tr","H0")

return(OUT)
}

