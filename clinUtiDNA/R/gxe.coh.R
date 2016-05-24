gxe.coh <- 
function(aff,unaff) {
cohort <- c(aff,unaff)
n <- sum(cohort)
X1=aff[1];  X2=aff[2]; X3=aff[3]; X4=aff[4];
Y1=unaff[1];  Y2=unaff[2]; Y3=unaff[3]; Y4=unaff[4];
int.est <- X1/(X1 + Y1) - X2/(X2 + Y2) - X3/(X3 + Y3) + X4/(X4 + Y4)
varint.est <- VarInt.coh(cohort)
lower <- int.est-qnorm(0.975)*varint.est^0.5
upper <- int.est+qnorm(0.975)*varint.est^0.5
res <- NULL
res$GxE.coh <- c(Estimate=int.est,Variance=varint.est,Lower=lower,Upper=upper)
return (res)
}
