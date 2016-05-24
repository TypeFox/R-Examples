PR2 <-
function(x,y,p,Rmat,rvec)
{

n=length(y);k=1
# Estimation of predictive regression and AR(p) model for the predictor
M = ARM(x,y,p)
b =M$b; a=M$a; cova=M$cova;covb=M$covb

if (min(Mod(polyroot(c(1,-a[1:p])))) <= 1) a[1:p] = ar.burg(x,aic=F,order.max=p)$ar

xmat0=M$xmat; ymat=M$ymat

# Bias-correction for AR(p) coefficients for k predictors
M = Shaman.Stine(x,p,a)
ac = M$coef;
vc = M$resid
mat = solve(M$mat)
covac = mat %*% cova %*% t(mat)

# Estimation of predictve regression, n: sample size, k; number of predictors
xmat <- cbind(xmat0,vc)
bc <- solve(crossprod(xmat)) %*% crossprod(xmat,ymat)
ec <- ymat - xmat %*% bc
s2 <- sum(ec^2)/(nrow(xmat)-length(bc))
covbc <- s2*solve(crossprod(xmat))

# Calculation of Covariance Matrix of AR(p) coefficients for k predictors
covcbc = bc[p+2]^2 * covac + covbc[2:(p+1),2:(p+1)]

# t-statistics and F-statistics
tstat = b/sqrt(diag(covb));  ptstat=2*(1-pt(abs(tstat),df=length(y)-length(b)))
if (p == 1) tstatc = bc[2:(p+1)]/sqrt(covcbc[1:p,1:p])
if (p > 1) tstatc = bc[2:(p+1)]/sqrt(diag(covcbc[1:p,1:p]))
tstatc = c(tstat[1],tstatc,bc[p+2]/sqrt(covbc[p+2,p+2]))
ptstatc=2*(1-pt(abs(tstatc),df=length(y)-length(bc)))

fstat= t(Rmat %*% matrix(b[2:(p+1),])-rvec) %*%   solve(Rmat %*% covb[2:(p+1),2:(p+1)] %*% t(Rmat))  %*% (Rmat%*%matrix(b[2:(p+1),])-rvec)/nrow(Rmat)
fstatc= t(Rmat %*% matrix(bc[2:(p+1),])-rvec) %*%   solve(Rmat %*% covcbc %*% t(Rmat))  %*% (Rmat%*%matrix(bc[2:(p+1),])-rvec)/nrow(Rmat)

fstatp=1-pf(fstat,df1=nrow(Rmat),df2=length(y)-nrow(b))
fstatcp= 1-pf(fstatc,df1=nrow(Rmat),df2=length(y)-nrow(bc))

#Names
LS=cbind(b,tstat,ptstat);IARM=cbind(bc,tstatc,ptstatc)

tem=character()
for(i in 1:k) {
tem1=paste("x",i,sep="")
for(j in 1:p) {
tem2=paste(tem1,-j,sep="("); tem3=paste(tem2,")",sep="")
tem=c(tem,tem3)
}}
tem=c("intercept",tem)
rownames(LS) <-tem;  colnames(LS) <- c("coefficients","t-stats","p-val")

tem4=paste("v",1:k,sep="")
rownames(IARM) <-c(tem,tem4);colnames(IARM) <- c("coefficients","t-stats","p-val")

tem1=c(paste("AR",1:p,sep=""),"intercept")
tem2=paste("x",1:k,sep="")

AR=a
rownames(AR)<- tem1
colnames(AR) <- tem2

AR.c=ac
rownames(AR.c)<- tem1
colnames(AR.c) <- tem2

Fstats=rbind(cbind(fstat,fstatc),cbind(fstatp,fstatcp))
rownames(Fstats)=c("F-stat","p-val"); colnames(Fstats)=c("OLS","Improved ARM")

return(list(LS=LS,IARM=IARM,AR=AR,ARc=AR.c,Fstats=Fstats,Covbc=covcbc))
}
