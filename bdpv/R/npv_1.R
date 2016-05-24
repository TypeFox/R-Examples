nPV <-
function(se, sp, prev, NPV0, PPV0, NPVpower=0.8, PPVpower=0.8, rangeP=c(0.05,0.95), nsteps=20, alpha=0.05, setnames=NULL)
{

arglist<-list(se=se, sp=sp, prev=prev, NPV0=NPV0, PPV0=PPV0, NPVpower=NPVpower, PPVpower=PPVpower)

arglength<-unlist(lapply(arglist, length))
maxlength<-max(arglength)

if(is.null(setnames)){setnames<-paste("Setting ", 1:maxlength, sep="")}
ARGLIST <- lapply(X=arglist, FUN=function(x){rep(x, length.out=maxlength)})
inDAT <- as.data.frame(ARGLIST)
rownames(inDAT) <- make.unique(rep(setnames, length.out=maxlength))
inPPV<-apply(X=inDAT, MARGIN=1, FUN=function(x){ppv(p=x["prev"], se=x["se"], sp=x["sp"])})
inNPV<-apply(X=inDAT, MARGIN=1, FUN=function(x){npv(p=x["prev"], se=x["se"], sp=x["sp"])})
CONF.LEVEL=1-alpha
Pseq <- seq(from=min(rangeP), to=max(rangeP), length.out=nsteps)

nlist<-list()

for (i in 1:maxlength)
{
PPV <- nPPV(propP=Pseq, se=ARGLIST$se[i], sp=ARGLIST$sp[i], prev=ARGLIST$prev[i],
PPV0=ARGLIST$PPV0[i], power=ARGLIST$PPVpower[i], conf.level=CONF.LEVEL )
NPV <- nNPV(propP=Pseq, se=ARGLIST$se[i], sp=ARGLIST$sp[i], prev=ARGLIST$prev[i],
NPV0=ARGLIST$NPV0[i], power=ARGLIST$NPVpower[i], conf.level=CONF.LEVEL )
nlist[[i]]<-list(NPV=NPV, PPV=PPV)
}

outDAT<-cbind(inDAT, trueNPV=inNPV, truePPV=inPPV)
RES<-list(outDAT=outDAT, nlist=nlist, NSETS=maxlength, nsteps=nsteps, rangeP=rangeP, propP=Pseq)
class(RES) <- "nPV"
return(RES)
}

