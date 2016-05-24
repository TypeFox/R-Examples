plotCI.UnlogCI<-function(x,...)
{

alt<-x$alternative

CIEs<-x$conf.int

rnames<-rownames(CIEs)

est<-CIEs[,"Estimate"]
lwr<-CIEs[,"lwr"]
upr<-CIEs[,"upr"]

names(est)<-rnames

args<-list(...)

args$estimate<-est
args$lower<-lwr
args$upper<-upr
args$alternative<-alt

do.call("plotCII", args)

}

