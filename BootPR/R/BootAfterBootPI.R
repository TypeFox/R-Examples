BootAfterBootPI <-
function(x,p,h,nboot,prob,type)
{x<-as.matrix(x)
if (type=="const")
M <- BootAfterBoot.PI(x,p,h,nboot,prob) 
if (type=="const+trend")
M <- BootAfterBootT.PI(x,p,h,nboot,prob) 
colnames(M$PI) <- paste(prob*100,"%",sep="");rownames(M$PI) <- paste("h",1:h,sep="")
colnames(M$forecast) <- "forecasts"; rownames(M$forecast) <- paste("h",1:h,sep="")
return(list(PI=M$PI,forecast=M$forecast))
}
