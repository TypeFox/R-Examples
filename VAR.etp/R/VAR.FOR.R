VAR.FOR <-
function(x,p,h,type="const",alpha=0.95)
{
M=VAR.est(x,p,type)
f=VAR.Fore(x,M$coef,p,h,type)
dims <- dim(f)
h <- dims[1]; k <- dims[2]
intv <- matrix(NA,h,2*k)
cr=qnorm(0.5*(1-alpha))

index=seq(from=1,by=2,length.out=k)
for(i in 1:h){
FMSE=VAR.mseh(x,p,h=i,type)
lower=f[i,]+cr*sqrt(diag(FMSE))
upper=f[i,]-cr*sqrt(diag(FMSE))
intv[i,index]=lower; intv[i,index+1]=upper
}

rownames(intv) <- rownames(f)
tem1 <- paste(colnames(f),"_lower",sep="")
tem2 <- paste(colnames(f),"_upper",sep="")
tem <- character(2*k)
tem[seq(1,(2*k-1),2)] <- tem1;tem[seq(2,(2*k),2)] <- tem2;
colnames(intv) <- tem
return(list(Intervals=intv,Forecast=f,Prob=alpha)) }
