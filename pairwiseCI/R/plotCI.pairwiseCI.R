plotCI.pairwiseCI<-function(x,...)
{

alt<-x$alternative

dat<-as.data.frame(x)

if(length(x$bynames)>1)
{
rnames<-paste(dat$by, dat$comparison, sep=": ")

}
else{
rnames<-dat$comparison
}

est<-dat$estimate
lwr<-dat$lower
upr<-dat$upper

names(est)<-names(lwr)<-names(upr)<-rnames

args<-list(...)

args$estimate<-est
args$lower<-lwr
args$upper<-upr
args$alternative<-alt

do.call("plotCII", args)

}




