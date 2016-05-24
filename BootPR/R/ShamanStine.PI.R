ShamanStine.PI <-
function(x,p,h,nboot,prob,type,pmax)
{
x<-as.matrix(x)
if (type=="const" & pmax==0)
M <- ShamanStine1.PI(x,p,h,nboot,prob)
if (type=="const+trend" & pmax==0)
M <- ShamanStine1T.PI(x,p,h,nboot,prob)
if (type=="const" & pmax > 0)
M <- ShamanStine2.PI(x,p,h,nboot,prob,pmax)
if (type=="const+trend" & pmax > 0)
M <- ShamanStine2T.PI(x,p,h,nboot,prob,pmax)
colnames(M$PI) <- paste(prob*100,"%",sep="");rownames(M$PI) <- paste("h",1:h,sep="")
colnames(M$forecast) <- "forecasts"; rownames(M$forecast) <- paste("h",1:h,sep="")
return(list(PI=M$PI,forecast=M$forecast))
}
