Andrews.Chen <-
function(x,p,h,type)
{
x<-as.matrix(x)
if (type=="const")
M <- Andrews.Chen1(x,p,h)
if (type=="const+trend")
M <- Andrews.Chen2(x,p,h)
rownames(M$coef) <- ARnames(p,type); colnames(M$coef) <- "coefficients"
rownames(M$ecmcoef) <- ARnames2(p,type); colnames(M$coef) <- "coefficients"
colnames(M$forecast) <- "forecasts"; rownames(M$forecast) <- paste("h",1:h,sep="")
return(list(coef=M$coef,ecm.coef=M$ecmcoef,resid=M$resid,forecast=M$forecast))
}
