bestfit.bootstrap<-function(best.lstsqr,boot,year,rawvar,rawcovar,s1,P,R,Q,covariate){

A<-best.lstsqr$A
B<-best.lstsqr$B
C<-best.lstsqr$C
mu_inf<-best.lstsqr$stationary.distribution$mean
V_inf<-best.lstsqr$stationary.distribution$covariance
E<-best.lstsqr$process.errors$residuals
sigma<-best.lstsqr$process.errors$sigma

# number of rounds to bootstrap
bootn<-boot

varind<-B!=0
covarind<-C!=0
if(length(covarind)<1) covarind<-NULL


# Define variates in which to place the bootstrapped values

bootB<-array(0,dim=c(P,P,bootn))
bootA<-matrix(0,bootn,P)
bootC<-array(0,dim=c(P,R,bootn))

bootsigma<-array(0,dim=c(P,P,bootn))

bootmuinf<-matrix(0,bootn,P)
bootVinf<-array(0,dim=c(P,P,bootn))

# Do the bootstrapping

for(jj in 1:bootn){

Estar<-t(E[sample(1:Q,replace=T),])


# reconstruct "raw" dataset using the covariates at time t
Xstar<-matrix(0,P,s1) # length of original data
kk<-1  # to step through the reconstructed residuals

for(ii in 1:s1){

if(ii==1) Xstar[,ii]<-t(rawvar[1,]) else {
if(year[ii]!=year[ii-1]) Xstar[,ii]<-t(rawvar[ii,])
if(year[ii]==year[ii-1]) {
	if(R>0){
	Xstar[,ii]<-A+B%*%Xstar[,ii-1]+C%*%as.numeric(rawcovar[ii,])+Estar[,kk]} else {
	Xstar[,ii]<-A+B%*%Xstar[,ii-1]+Estar[,kk]}
	kk<-kk+1}
}

}

Xstar<-t(Xstar)
XX<-cbind(Xstar[1:s1-1,],Xstar[2:s1,])

# throw out non-overalpping years
XX<-XX[which(year[1:s1-1]==year[2:s1]),]
lagstateS<-XX[,1:P,drop=F]
statevarS<-XX[,-c(1:P),drop=F]

# keep the save covariate matrices as used in the main fits

# initialize variates for fitting the bootstrapped model
Es<-matrix(0,Q,P)
Astar<-matrix(0,P,1)
Bstar<-matrix(0,P,P)
Cstar<-matrix(0,P,R)

# loop over each variate, fitting least squares estimates as we go
for(dv in 1:P){
Ys<-statevarS[,dv]
lv<-which(varind[dv,]==1)
cv<-which(covarind[dv,]==1)
nvar<-sum(varind[dv,])
ncovar<-sum(covarind[dv,])

Xs<-cbind(1,lagstateS[,lv],covariate[,cv])

# least squares estimates
betastar<-solve( (t(Xs)%*%Xs), (t(Xs)%*%Ys) )

# calculate the residuals for that dependent variate
Es[,dv]<-Ys-Xs%*%betastar

# parameter estimates
Astar[dv,]<-betastar[1,1]
if (length(lv)>0) Bstar[dv,lv]<-t(betastar[2:(nvar+1),1])
if (length(cv)>0) Cstar[dv,cv]<-t(betastar[(nvar+2):length(betastar[,1]),1])
}

# CALCULATE & SAVE WHAT'S BEING BOOTSTRAPPED
bootA[jj,]<-t(Astar)
bootB[1:P,1:P,jj]<-Bstar
bootC[1:P,1:R,jj]<-Cstar

# variance-covariance matrix of the error terms
sigma1<-t(Es)%*%E/Q
bootsigma[1:P,1:P,jj]<-sigma1

# mu_inf (= mean of that stationary distribution)
bootmuinf[jj,]<-solve(diag(P)-Bstar)%*%Astar

# V_inf (= variance-covariance of the stationary distribution)
vecV<-solve(diag(P*P)-kronecker(Bstar,Bstar))%*%as.vector(sigma1)
bootVinf[1:P,1:P,jj]<-matrix(vecV,P,P)

}

# BOOTSTRAP OUTPUT

# get upper and lower values for A
lA<-matrix(0,P,1); uA<-matrix(0,P,1)
for(i in 1:P){
	sorted<-sort(bootA[,i])
	lA[i,1]<-sorted[floor(0.025*bootn)]
	uA[i,1]<-sorted[bootn-floor(0.025*bootn)]}

# get upper and lower values for the B matrix
lB<-matrix(0,P,P); uB<-matrix(0,P,P)
for(i in 1:P){
for(j in 1:P){
	sorted<-sort(bootB[i,j,])
	lB[i,j]<-sorted[floor(0.025*bootn)]
	uB[i,j]<-sorted[bootn-floor(0.025*bootn)]}}

# get upper and lower values for the C matrix
if(R>0){
lC<-matrix(0,P,R); uC<-matrix(0,P,R)
for(i in 1:P){
for(j in 1:R){
	sorted<-sort(bootC[i,j,])
	lC[i,j]<-sorted[floor(0.025*bootn)]
	uC[i,j]<-sorted[bootn-floor(0.025*bootn)]}}
} else {lC<-NULL; uC<-NULL}

# get upper and lower values for the sigma matrix
lSIGMA<-matrix(0,P,P); uSIGMA<-matrix(0,P,P)
for(i in 1:P){
for(j in 1:P){
	sorted<-sort(bootsigma[i,j,])
	lSIGMA[i,j]<-sorted[floor(0.025*bootn)]
	uSIGMA[i,j]<-sorted[bootn-floor(0.025*bootn)]}}

# get upper and lower values for the Vinf matrix
lVinf<-matrix(0,P,P); uVinf<-matrix(0,P,P)
for(i in 1:P){
for(j in 1:P){
	sorted<-sort(bootVinf[i,j,])
	lVinf[i,j]<-sorted[floor(0.025*bootn)]
	uVinf[i,j]<-sorted[bootn-floor(0.025*bootn)]}}

# get the upper and lower values for the mean of the stationary distribution
lmu_inf<-matrix(0,P,1); umu_inf<-matrix(0,P,1)
for(i in 1:P){
	sorted<-sort(bootmuinf[,i])
	lmu_inf[i,1]<-sorted[floor(0.025*bootn)]
	umu_inf[i,1]<-sorted[bootn-floor(0.025*bootn)]}

boot.A<-A
boot.A.ind<-which(lA<0&uA>0)
boot.A[boot.A.ind]<-0

boot.B<-B
boot.B.ind<-which(lB<0&uB>0)
boot.B[boot.B.ind]<-0

boot.C<-C
boot.C.ind<-which(lC<0&uC>0)
boot.C[boot.C.ind]<-0

boot.sigma<-sigma
boot.sigma.ind<-which(lSIGMA<0&uSIGMA>0)
boot.sigma[boot.sigma.ind]<-0

boot.mu_inf<-mu_inf
boot.mu_inf.ind<-which(lmu_inf<0&umu_inf>0)
boot.mu_inf[boot.mu_inf.ind]<-0

boot.V_inf<-V_inf
boot.V_inf.ind<-which(lVinf<0&uVinf>0)
boot.V_inf[boot.V_inf.ind]<-0


list(
	lowerA	=	lA,
	upperA	=	uA,
	bootA	=	boot.A,
	lowerB	=	lB,
	upperB	=	uB,
	bootB	=	boot.B,
	lowerC	=	lC,
	upperC	=	uC,
	bootC	=	boot.C,
	stationary.distribution=list(
		lower.mean	=	lmu_inf,
		upper.mean	=	umu_inf,
		boot.mean	=	boot.mu_inf,
		lower.covariance	=	lVinf,
		upper.covariance	=	uVinf,
		boot.covariance		=	boot.V_inf	),
	process.errors=list(
		lower.sigma	=	lSIGMA,
		upper.sigma	=	uSIGMA,
		boot.sigma	=	boot.sigma	)
	)

}


