bestfit.lstsqr<-function(statevar,lagstate,covariate,P,R,Q,bestGlobalB,bestGlobalC,
								namesvar,namescovar){
#====================================================================================
# RECALCULATE THE LEAST-SQUARES ESTIMATES FOR THE SELECTED BEST-FIT MODEL:
#====================================================================================

# Reset interaction restrictions to reconstruct the best-fit model
varind<-bestGlobalB!=0
covarind<-bestGlobalC!=0
if(length(covarind)<1) covarind<-NULL


# initialize variates
E<-matrix(0,nrow=Q,ncol=P)	# residuals
A<-matrix(0,nrow=P,ncol=1)	# intercepts for the variates
B<-matrix(0,nrow=P,ncol=P)	# parameters for the variates
C<-matrix(0,nrow=P,ncol=R)	# parameters for the covariates
Yhat<-matrix(0,nrow=Q,ncol=P)	# estimates
varY<-matrix(0,nrow=Q,ncol=P)	# variates for total variance
varDY<-matrix(0,nrow=Q,ncol=P)	# variates for per capita variance


for(dv in 1:P){		#  ************* loop over each variate *************

Y<-statevar[,dv]
lv<-which(varind[dv,]==1)
cv<-which(covarind[dv,]==1)

nvar<-sum(varind[dv,])

X<-cbind(1,lagstate[,lv],covariate[,cv])		# predictors

# least squares estimates
beta<-solve( (t(X)%*%X),(t(X)%*%Y) )

# calculate residuals for the dependent variate
Yhat[,dv]<-X%*%beta
E[,dv]<-Y-Yhat[,dv]

# parameter estimates
A[dv,]<-beta[1,1]
if (length(lv)>0) B[dv,lv]<-beta[2:(nvar+1),1]
if (length(cv)>0) C[dv,cv]<-beta[(nvar+2):length(beta[,1]),1]

# accumulate Y variates for calculation of R^2 for observed vs. predicted
varY[,dv]<-Y-mean(Y)

# accumulate Y variates for calculation of R^2 for observed change in x vs. predicted
varDY[,dv]<-Y-lagstate[,dv]

}				#  ************* END loop *************


# ASSESS THE LEAST-SQUARES FITS

# calculate LS log-likelihood
sigma<-t(E)%*%E/Q
lnlike<- -Q*(P/2)*log(2*pi)-(Q/2)*log(det(sigma))-Q*P/2

# calculate LS explained variance
varMatrix<-t(varY)%*%varY/Q
R2<-1-diag(sigma)/diag(varMatrix)

# calculate LS explained variance in change in x
varMatrix_D<-t(varDY)%*%varDY/Q
R2_D<-1-diag(sigma)/diag(varMatrix_D)

# correlation matrix for noise
d<-diag(1/sqrt(diag(sigma)))
corrmatrix<-d%*%sigma%*%d

# parameter count = number of non-zero parameters in A, B, C, and some of sigma
par<-sum(A!=0,B!=0,C!=0)+P*(P+1)/2

# AIC and BIC
lsAIC<- -2*lnlike+2*par
lsBIC<- -2*lnlike+par*log(Q)

# CALCULATE THE STATIONARY DISTRIBUTION

mu_inf<-solve(diag(P)-B)%*%A		# mean of stationary distribution
vecV<-solve(diag(P*P)-kronecker(B,B))%*%as.vector(sigma)
V_inf<-matrix(vecV,nrow=P,ncol=P)	# variance-covariance of stationary distribution


# RESULTS

rownames(A)<-namesvar;colnames(A)<-"a"
dimnames(B)<-list(namesvar,namesvar)
dimnames(C)<-list(namesvar,namescovar)
R2.values<-cbind(R2,R2_D)
rownames(R2.values)<-namesvar
dimnames(mu_inf)<-list(namesvar,"mu")
dimnames(V_inf)<-list(namesvar,namesvar)
colnames(E)<-namesvar
dimnames(sigma)<-list(namesvar,namesvar)
dimnames(corrmatrix)<-list(namesvar,namesvar)

list(
A		=	A,
B		=	B,
C		=	C,
log.likelihood		=	lnlike,
AIC					=	lsAIC,
BIC					=	lsBIC,
R2.values			=	R2.values,
stationary.distribution=list(
		mean		=	mu_inf,
		covariance	=	V_inf		),
process.errors=list(
		residuals	=	E,
		covariance	=	sigma,
		corrmatrix	=	corrmatrix	)
)

}
