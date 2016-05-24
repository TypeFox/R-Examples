Spec.shape <-
function(x)
{
iad <- numeric(11)
icvm <- numeric(11)
ima <- numeric(11)

for (j in 1:11)
{
iacm <- IACM( (j-1)*0.1,x)
iad[j] <- iacm[1]
icvm[j] <- iacm[2]
ima[j] <- iacm[3]
}

ad <- ISIMP(0,1,iad)
cvm <- ISIMP(0,1,icvm)
ma <- ISIMP(0,1,ima)
return(list(AD=ad,CVM=cvm,M=ma))
}
