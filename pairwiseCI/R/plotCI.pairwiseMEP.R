`plotCI.pairwiseMEP` <-
function(x, whichep=NULL, ...)
{

pargs<-list(...)

tabm<-table(x$method)

if(is.null(whichep))
{
# plot all in one plot!

if(length(tabm)>1)
{
warning("Plotting all intervals in the same plot makes hardly sense, when different methods have been applied for calculation!")}
}

dat<-as.data.frame(x, whichep=whichep)

EST<-dat$estimate
LOW<-dat$lower
UPP<-dat$upper

if("by" %in% names(dat))
{
names(EST)<-paste(dat$response, dat$by, dat$comparison, sep=": ")
}
else{
names(EST)<-paste(dat$response, dat$comparison, sep=": ")
}

pargs$estimate<-EST
pargs$lower<-LOW
pargs$upper<-UPP
pargs$alternative<-x$alternative

do.call("plotCII", pargs)

}

