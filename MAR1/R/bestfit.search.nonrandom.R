bestfit.search.fwdstep<-function(statevar,lagstate,covariate,P,R,Q,indexBGlobal,indexCGlobal,...){
#====================================================================================
# INITIALIZE BEST-FIT MODEL SEARCH (lowest AIC):
#====================================================================================

bestGlobalB<-matrix(0,nrow=P,ncol=P)
bestGlobalC<-matrix(0,nrow=P,ncol=R)

nmodels<-0

for(dv in 1:P){	#  *************** START Main loop ***************
# (for each of the variates...)

indexB<-indexBGlobal[dv,]
indexC<-indexCGlobal[dv,]
indexBC<-c(indexB,indexC)

bestAIC	<-10000000
lowestAIC	<- 1000000

while(lowestAIC<bestAIC){	# ^^^^^^ START adding loop ^^^^^^
# (find bestAIC by adding variables 1 by 1 to lowestAIC models)

bestAIC<-lowestAIC

stilltestBC<-which(indexBC==0.5)
testing.orig<-rep(0,length(stilltestBC))

for(rept in 0:length(stilltestBC)){			# ------ START replacement loop ------
# (find lowestAIC of models built by cycling '1' through 0.5 index values)

E<-matrix(0,nrow=Q,ncol=1)	# residuals
A<-matrix(0,nrow=1,ncol=1)	# intercepts for the variates
B<-matrix(0,nrow=1,ncol=P)	# parameters for the variates
C<-matrix(0,nrow=1,ncol=R)	# parameters for the covariates
Yhat<-matrix(0,nrow=Q,ncol=1)	# estimates

testing<-testing.orig
testing[rept]<-1

varcovarind<-indexBC
varcovarind[varcovarind==0.5]<-testing
# varind = matrix for the variates : which interactions to include
varind<-varcovarind[1:P]
# covarind = matrix for the covariates : which interactions to include
if(R>0) covarind<-varcovarind[(P+1):length(varcovarind)] else covarind<-NULL

Y<-statevar[,dv]						# dependent variate
lv<-which(varind==1)
cv<-which(covarind==1)

nvar<-sum(varind)

X<-cbind(1,lagstate[,lv],covariate[,cv])		# predictors

# least squares estimates
beta<-solve( (t(X)%*%X),(t(X)%*%Y) )

# calculate residuals for the dependent variate
Yhat<-X%*%beta
E[,1]<-Y-Yhat

# parameter estimates
A<-beta[1,1]
if (length(lv)>0) B[1,lv]<-beta[2:(nvar+1),1]
if (length(cv)>0) C[1,cv]<-beta[(nvar+2):length(beta[,1]),1]

# calculate log-likelihood
sigma<-t(E)%*%E/Q
lnlike<- -Q*(1/2)*log(2*pi)-(Q/2)*log(det(sigma))-Q/2

par<-sum(A!=0,B!=0,C!=0)+1
SSAIC<- -2*lnlike+2*par

if(SSAIC<lowestAIC){
	lowestAIC<-SSAIC
	lowestB<-B
	lowestC<-C
	}

nmodels<-nmodels+1

}						# ------- END replacement loop -------

if(lowestAIC<bestAIC){
	bestB<-lowestB
	bestC<-lowestC
	
	indexB[which(bestB!=0)]<-1
	indexC[which(bestC!=0)]<-1
	indexBC<-c(indexB,indexC)
	}


}					# ^^^^^^^ END adding loop ^^^^^^^

bestGlobalB[dv,]<-bestB
bestGlobalC[dv,]<-bestC

}			#  **************** END Main loop ****************


list(bestGlobalB=bestGlobalB,bestGlobalC=bestGlobalC,all.models=NULL)

}

