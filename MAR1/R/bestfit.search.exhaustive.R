bestfit.search.exhaustive<-function(statevar,lagstate,covariate,P,R,Q,indexBGlobal,indexCGlobal,ntop){
#====================================================================================
# INITIALIZE BEST-FIT MODEL SEARCH (lowest AIC):
#====================================================================================

#library(leaps)

bestGlobalB<-matrix(0,nrow=P,ncol=P)
bestGlobalC<-matrix(0,nrow=P,ncol=R)

if(ntop==F) ntop<-2
all.models<-vector("list",P)

for(dv in 1:P){	#  *************** START Main loop ***************
# (for each of the variates...)

indexB<-indexBGlobal[dv,]
indexC<-indexCGlobal[dv,]
indexBC<-c(indexB,indexC)

E<-matrix(0,nrow=Q,ncol=1)	# residuals
A<-matrix(0,nrow=1,ncol=1)	# intercepts for the variates
B<-matrix(0,nrow=1,ncol=P)	# parameters for the variates
C<-matrix(0,nrow=1,ncol=R)	# parameters for the covariates
Yhat<-matrix(0,nrow=Q,ncol=1)	# estimates

# varind = matrix for the variates : which interactions to include
varind<-indexB
varind[varind==.5]<-1

# covarind = matrix for the covariates : which interactions to include
if(R>0) {
covarind<-indexC
covarind[covarind==.5]<-1} else {
covarind<-NULL}

Y<-statevar[,dv]						# dependent variate
lv<-which(varind==1)
cv<-which(covarind==1)

nvar<-sum(varind)

X<-cbind(1,lagstate[,lv],covariate[,cv])		# predictors
X<-as.data.frame(X[,-1])

best.subsets<-summary(regsubsets(x=X,y=Y,nbest=round(ntop/2),nvmax=ncol(X),
	force.in=which(indexBC[indexBC!=0]==1)))$which[,-1]
best.subsets<-best.subsets[,order(as.numeric(substr(colnames(best.subsets),2,5)))]
if(any(indexBC==1)){
best.subsets<-rbind(NA,best.subsets)
best.subsets[1,]<-indexBC[indexBC!=0]==1	}

rownames(best.subsets)<-NULL
sub.aics<-data.frame(AIC=TRUE,best.subsets)

for(i in 1:nrow(best.subsets)){
Xt<-X[,best.subsets[i,]]
sub.glm<-glm(Y~.,data=data.frame(Xt))
sub.aics[i,which(sub.aics[i,]==T)]<-c(AIC(sub.glm),sub.glm$coefficients[-1])
}

sub.aics<-sub.aics[order(sub.aics$AIC),]
best.subset<-sub.aics[1,-1]
indexBC[indexBC!=0]<-as.numeric(best.subset)
indexB<-indexBC[1:P]
indexC<-indexBC[(P+1):length(indexBC)]

bestGlobalB[dv,]<-as.numeric(indexB)
bestGlobalC[dv,]<-as.numeric(indexC)

all.models[[dv]]<-matrix(c(999,indexBGlobal[dv,],indexCGlobal[dv,]),
	nrow=nrow(sub.aics),ncol=P+R+1,byrow=T)
all.models[[dv]][all.models[[dv]]!=0]<-as.matrix(sub.aics)

}			#  **************** END Main loop ****************


list(bestGlobalB=bestGlobalB,bestGlobalC=bestGlobalC,all.models=all.models)

}

