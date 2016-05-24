MandelPvalue <-
function(hfobj)
{
a <- nlevels(hfobj$tall$rows)
b <- nlevels(hfobj$tall$cols)
ymtx <- matrix(hfobj$tall$y,nrow=a,ncol=b,byrow=T)
coldevs <- apply(ymtx,2,mean)-mean(ymtx)
rowdevs <- apply(ymtx,1,mean)-mean(ymtx)
SSRow <- b*sum(rowdevs^2)
SSCol <- a*sum(coldevs^2)
SSTot <- (a*b-1)*var(hfobj$tall$y)
slopes <- ymtx %*% coldevs/sum(coldevs^2)
# SSMandel <- (t(slopes) %*% rowdevs)^2/sum(rowdevs^2) * sum(coldevs^2)
SSMandel <- sum((slopes-1)^2) * sum(coldevs^2)
SSE <- SSTot-SSMandel-SSRow-SSCol
dfE <- ((a-1)*(b-2))
MSE <- SSE/dfE
Fratio <- (SSMandel/(a-1))/MSE
pvalue <- 1-pf(Fratio,(a-1),dfE)
SumSq <- c(SSRow=SSRow,SSCol=SSCol,SSMandel=SSMandel,SSE=SSE,SSTot=SSTot)
list(pvalue=pvalue,SumSq=SumSq,Fratio=Fratio,df=c(a-1,dfE))
#list(mandel.pvalue=pvalue)
}
