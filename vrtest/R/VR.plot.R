VR.plot <-
function(y,kvec)
{
val <- matrix(NA,nrow=max(kvec),ncol=3)
for( i in 2:max(kvec))
{
tem1 <- stat.plot(y,i)$vr
tem2 <- stat.plot(y,i)$se
val[i,] <- c(tem1,1-1.96*tem2,1+1.96*tem2)}

matplot(val,type="l",col=c(2,4,4),xlab="holding period",ylab="variance ratio",lwd=c(5,2,2))
abline(h=1)
grid(nx=max(kvec),lwd=1)
title(main = "Variance Ratios and 95% confidence band")
VAL <- as.matrix(val[2:max(kvec),1])
rownames(VAL) <- paste("k=",2:max(kvec),sep="")
colnames(VAL) <- "VR"
return(list(VR=VAL))
}
