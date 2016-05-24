LS.AR <-
function(x,p,h,type,prob)
{
x<-as.matrix(x)
if (type=="const")
M <- OLS.AR(x,p,h,prob)
if (type=="const+trend")
M <- OLS.ART(x,p,h,prob)
rownames(M$coef) <- ARnames(p,type); colnames(M$coef) <- "coefficients"
colnames(M$PI) <- paste(prob*100,"%",sep="");rownames(M$PI) <- paste("h",1:h,sep="")
colnames(M$forecast) <- "forecasts"; rownames(M$forecast) <- paste("h",1:h,sep="")
return(list(coef=M$coef,resid=M$resid,forecast=M$forecast,PI=M$PI))
}
