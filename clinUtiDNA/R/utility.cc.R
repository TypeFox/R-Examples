utility.cc <- 
function(cases, contr, pD,pG) {
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
if(pG[1]<1) {
  uti.est <- pG*(1-pG)*int.est
  varuti.est <- varint.est
}
if(pG[1]>1) {
  XG=pG[1]; YG=pG[2]
  nG=XG+YG
  uti.est <- XG/nG*(1-XG/nG)*int.est
  if(pD[1]<1) varuti.est <- VarUti.cc.pD(pD,c(cases,contr,pG))
  if(pD[1]>1) varuti.est <- VarUti.cc(c(cases,contr,pD,pG))
}
lower <- uti.est-qnorm(0.975)*varuti.est^0.5
upper <- uti.est+qnorm(0.975)*varuti.est^0.5
res <- NULL
res$Utility.cc <- c(Estimate=uti.est,Variance=varuti.est,Lower=lower,Upper=upper)
return (res)
}
