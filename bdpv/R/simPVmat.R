simPVmat <-
function(se, sp, pr, n1, n0, NPV0, PPV0, conf.level=0.95, NSIM=500, setnames=NULL)
{

DAT<-cbind(se, sp, pr, NPV0, PPV0, n1, n0)
colnames(DAT)<-c("se", "sp", "pr", "NPV0", "PPV0", "n1", "n0")

inDAT <- as.data.frame(DAT)

if(is.null(setnames))
{SNAMES<-paste("Setting ", 1:nrow(inDAT), sep="")}
else{SNAMES<-rep(setnames, length.out=nrow(inDAT))}

OUTNPV<-NULL
OUTPPV<-NULL

for(i in 1:nrow(inDAT))
{

OUT<-simPV(se=inDAT[i,"se"], sp=inDAT[i,"sp"], pr=inDAT[i,"pr"],
 n1=inDAT[i,"n1"], n0=inDAT[i,"n0"], NPV0=inDAT[i,"NPV0"], PPV0=inDAT[i,"PPV0"],
 conf.level=conf.level, NSIM=NSIM)


OUTNPV<-rbind(OUTNPV, OUT["NPV",])
OUTPPV<-rbind(OUTPPV, OUT["PPV",])

}

inDAT$n <- inDAT$n0+inDAT$n1

rownames(inDAT)<-rownames(OUTNPV)<-rownames(OUTPPV)<-SNAMES

return(list(INDAT=inDAT, NPV=OUTNPV, PPV=OUTPPV, NSIM=NSIM, conf.level=conf.level))
}

