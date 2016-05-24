PR3 <-
function(x,y,p,Rmat,rvec)
{

# Estimation of predictve regression, n: sample size, k; number of predictors
n=length(y); k = ncol(x)
ymat <- matrix(y[(p+1):n])
xmat0 <- matrix(1,nrow=n-p,k*p+1)
index1 = p:(n-1); 
if (p==1) index2=2:(k+1)
if (p>1) index2=seq.int(2,k*p,p)
for (i in 1:p){
xmat0[,index2] <- x[index1,];
index1=index1-1; index2=index2+1
}
b <- solve(crossprod(xmat0)) %*% crossprod(xmat0,ymat)
e <- ymat - xmat0 %*% b
s2 <- sum(e^2)/(nrow(xmat0)-length(b))
covb <- s2*solve(crossprod(xmat0))

# Estimation of k predctors, AR(p) model
amat = matrix(NA,nrow=p+1,ncol=k)
covamat = matrix(NA,nrow=p,ncol=p*k)
emat=matrix(NA,nrow=n-p,k)
x_mat = matrix(NA,nrow=p+1,ncol=k)
index = 1:p; 
for (i in 1:k){
A <- OLS.AR(x[,i],p)
amat[,i] <- A$coef; covamat[,index] <- A$covmat[1:p,1:p]; emat[,i] <- A$resid
tem1=A$xmat; tem2=A$resid
x_mat[,i]=solve(t(tem1) %*% tem1) %*% t(tem1) %*% tem2
index=index+p; 
}


# Bias-correction for AR(p) coefficients for k predictors
acmat = matrix(NA,nrow=p+1,ncol=k)
covacmat = matrix(NA,nrow=p,ncol=p*k)
vcmat=matrix(NA,nrow=n-p,k)
index = 1:p
for (i in 1:k){

if (min(Mod(polyroot(c(1,-amat[1:p,i])))) <= 1) amat[1:p,i] = ar.burg(x[,i],aic=F,order.max=p)$ar

M = Shaman.Stine(x[,i],p,amat[,i])
acmat[,i] = M$coef;
vcmat[,i] = M$resid
mat = solve(M$mat)
covacmat[,index] = mat %*% covamat[,index] %*% t(mat)
index=index+p
}

# Calculation of Covariance Matrix of AR(p) coefficients for k predictors
covac_mat=matrix(NA,nrow=p*k,ncol=p*k)
index1=1:p; index=1:p
for( i in 1:k)
{index2=1:p
for( j in 1:k){
if(i==j) {covac_mat[index1,index2] = covacmat[,index]; index=index+p; index2=index2+p }
if(i != j) {
tem1=x_mat[1:p,i];tem2=x_mat[1:p,j]; 
covac_mat[index1,index2] = mat %*% as.matrix(tem1,nrow=p) %*% matrix(tem2,nrow=1) %*% t(mat);
index2=index2+p}
}
index1=index1+p
}

# Augmented Regression for the Predtive Regression
xmat <- cbind(xmat0,vcmat)
bc <- solve(crossprod(xmat)) %*% crossprod(xmat,ymat)
ec <- ymat - xmat %*% bc
s2 <- sum(ec^2)/(nrow(xmat)-length(bc))
covbc <- s2*solve(crossprod(xmat))
bcphi= rev(rev(bc)[1:k])

tem1=matrix(NA,nrow=k,ncol=k)
for(i in 1:k){
for(j in 1:k){
tem1[i,j]=bcphi[i]*bcphi[j]
}}
tem2 = tem1 %x% matrix(1,nrow=p,ncol=p)


# Bias-reduced Covariance matrix for the coefficients in Predictive Regression
covcbc=tem2*covac_mat + covbc[2:(p*k+1),2:(p*k+1)]

# t-statistics
tstat = b/sqrt(diag(covb));  ptstat=2*(1-pt(abs(tstat),df=length(y)-length(b)))
tstatc = bc[2:(p*k+1)]/sqrt(diag(covcbc))
tstatc = c(tstat[1],tstatc,bc[(k*p+2):(k*p+1+k)]/sqrt(diag(covbc[(k*p+2):(k*p+1+k),(k*p+2):(k*p+1+k)])))
ptstatc=2*(1-pt(abs(tstatc),df=length(y)-length(bc)))

# F-statistics
fstat= t(Rmat %*% matrix(b[2:(p*k+1),])-rvec) %*%   solve(Rmat %*% covb[2:(p*k+1),2:(p*k+1)] %*% t(Rmat))  %*% (Rmat%*%matrix(b[2:(p*k+1),])-rvec)/nrow(Rmat)
fstatc= t(Rmat %*% matrix(bc[2:(p*k+1),])-rvec) %*%   solve(Rmat %*% covcbc %*% t(Rmat))  %*% (Rmat%*%matrix(bc[2:(p*k+1),])-rvec)/nrow(Rmat)

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
rownames(LS) <-tem; colnames(LS) <- c("coefficients","t-stats","p-val")

tem4=paste("v",1:k,sep="")
rownames(IARM) <-c(tem,tem4); colnames(IARM) <- c("coefficients","t-stats","p-val")

tem1=c(paste("AR",1:p,sep=""),"intercept")
tem2=paste("x",1:k,sep="")

AR=amat
rownames(AR)<- tem1
colnames(AR) <- tem2

AR.c=acmat
rownames(AR.c)<- tem1
colnames(AR.c) <- tem2

Fstats=rbind(cbind(fstat,fstatc),cbind(fstatp,fstatcp))
rownames(Fstats)=c("F-stat","p-val"); colnames(Fstats)=c("OLS","Improved ARM")

return(list(LS=LS,IARM=IARM,AR=AR,ARc=AR.c,Fstats=Fstats,Covbc=covcbc))
}
