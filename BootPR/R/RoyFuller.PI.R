RoyFuller.PI <-
function(x,p,h,nboot,prob,type,pmax)
{
x<-as.matrix(x)
if (type=="const" & pmax==0)
M <- RoyFuller1.PI(x,p,h,nboot,prob,"const")
if (type=="const+trend" & pmax==0)
M <- RoyFuller1T.PI(x,p,h,nboot,prob,"const+trend")
if (type=="const" & pmax > 0)
M <- RoyFuller2.PI(x,p,h,nboot,prob,"const",pmax)
if (type=="const+trend" & pmax > 0)
M <- RoyFuller2T.PI(x,p,h,nboot,prob,"const+trend",pmax)
colnames(M$PI) <- paste(prob*100,"%",sep="");rownames(M$PI) <- paste("h",1:h,sep="")
colnames(M$forecast) <- "forecasts"; rownames(M$forecast) <- paste("h",1:h,sep="")
return(list(PI=M$PI,forecast=M$forecast))
}
