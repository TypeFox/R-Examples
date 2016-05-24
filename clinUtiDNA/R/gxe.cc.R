gxe.cc <- 
function(cases, contr, pD) {
n1 <- sum(cases)
n0 <- sum(contr)
X1=cases[1];  X2=cases[2]; X3=cases[3]; X4=cases[4];
Y1=contr[1];  Y2=contr[2]; Y3=contr[3]; Y4=contr[4];
res <- NULL
if(pD[1]<1) {
q <- pD*n0/(1-pD)/n1
int.est <- q*(X1/(q*X1 + Y1) - X2/(q*X2 + Y2) - X3/(q*X3 + Y3) + X4/(q*X4 + Y4))
varint.est <- VarInt.cc.pD(pD,c(cases,contr))
}
if(pD[1]>1) {
XD=pD[1];YD=pD[2]
nD=XD+YD
q <- XD/nD*n0/(YD/nD)/n1
int.est <-  q*(X1/(q*X1 + Y1) - X2/(q*X2 + Y2) - X3/(q*X3 + Y3) + X4/(q*X4 + Y4))
varint.est <- VarInt.cc(c(cases,contr,pD))
}
lower <- int.est-qnorm(0.975)*varint.est^0.5
upper <- int.est+qnorm(0.975)*varint.est^0.5
res <- NULL
res$GxE.cc <- c(Estimate=int.est,Variance=varint.est,Lower=lower,Upper=upper)
return (res)
}
