top.bestfit<-function(ntop,all.models,
	statevar,lagstate,covariate,P,R,Q,namesvar,namescovar){

sortAIC<-function(x){x<-x[order(x[,1]),]}

u.models<-lapply(all.models,unique)
u.models<-lapply(u.models,sortAIC)
u.models<-lapply(u.models,matrix,ncol=ncol(all.models[[1]]))

if(ntop>sum(unlist(lapply(u.models,nrow)))) ntop<-sum(unlist(lapply(u.models,nrow)))

top.models<-array(NA,dim=c(length(u.models),ncol(u.models[[1]])-1,ntop))
dimnames(top.models)<-list(namesvar,c(namesvar,namescovar),1:ntop)

for(i in 1:nrow(top.models)) top.models[i,,1]<-u.models[[i]][1,-1]

next.best<-u.models
for(i in 1:length(u.models)) {
next.best[[i]]<-matrix(u.models[[i]][-1,-1],ncol=ncol(top.models))
}

################################################################
# CALCULATE THE LEAST-SQUARES ESTIMATES FOR THE BEST-FIT MODEL:

tmp.model<-top.models[,,1,drop=F]
dim(tmp.model)<-dim(tmp.model)[1:2]

varind<-tmp.model[,1:nrow(tmp.model),drop=F]!=0
if(ncol(tmp.model)==nrow(tmp.model)) covarind<-NULL else {
covarind<-tmp.model[,-c(1:nrow(tmp.model)),drop=F]!=0 }

# initialize variates
E<-matrix(0,nrow=Q,ncol=P)	# residuals
A<-matrix(0,nrow=P,ncol=1)	# intercepts for the variates
B<-matrix(0,nrow=P,ncol=P)	# parameters for the variates
C<-matrix(0,nrow=P,ncol=R)	# parameters for the covariates
Yhat<-matrix(0,nrow=Q,ncol=P)	# estimates

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

}				#  ************* END variate loop *************

# ASSESS THE LEAST-SQUARES FITS

# calculate LS log-likelihood
sigma<-t(E)%*%E/Q
lnlike<- -Q*(P/2)*log(2*pi)-(Q/2)*log(det(sigma))-Q*P/2

# parameter count = number of non-zero parameters in A, B, C, and some of sigma
par<-sum(A!=0,B!=0,C!=0)+P*(P+1)/2

dimnames(top.models)[[3]][1]<-(-2*lnlike+2*par)
################################################################

z<-1
if(sum(unlist(lapply(next.best,nrow)))>0){
for(z in 2:ntop){		#  ========= START n top model loop =========

AICS<-NULL

for(i in 1:nrow(top.models)){		#  ^^^^^^^^^^ START row replacement loop ^^^^^^^^^^

if(nrow(next.best[[i]])==0){
	AICS<-c(AICS,NA)
	next}

tmp.model<-top.models[,,1,drop=F]
dim(tmp.model)<-dim(tmp.model)[1:2]
tmp.model[i,]<-next.best[[i]][1,]

varind<-tmp.model[,1:nrow(tmp.model),drop=F]!=0
if(ncol(tmp.model)==nrow(tmp.model)) covarind<-NULL else {
covarind<-tmp.model[,-c(1:nrow(tmp.model)),drop=F]!=0 }

################################################################
# CALCULATE THE LEAST-SQUARES ESTIMATES FOR THE TEMPORARY MODEL:

# initialize variates
E<-matrix(0,nrow=Q,ncol=P)	# residuals
A<-matrix(0,nrow=P,ncol=1)	# intercepts for the variates
B<-matrix(0,nrow=P,ncol=P)	# parameters for the variates
C<-matrix(0,nrow=P,ncol=R)	# parameters for the covariates
Yhat<-matrix(0,nrow=Q,ncol=P)	# estimates

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

}				#  ************* END variate loop *************

# ASSESS THE LEAST-SQUARES FITS

# calculate LS log-likelihood
sigma<-t(E)%*%E/Q
lnlike<- -Q*(P/2)*log(2*pi)-(Q/2)*log(det(sigma))-Q*P/2

# parameter count = number of non-zero parameters in A, B, C, and some of sigma
par<-sum(A!=0,B!=0,C!=0)+P*(P+1)/2

AICS<-c(AICS,(-2*lnlike+2*par))
################################################################

}						#  ^^^^^^^^^^ END row replacement loop ^^^^^^^^^^

top.models[,,z]<-top.models[,,1]
minAIC<-which(AICS==min(AICS,na.rm=T))
top.models[minAIC,,z]<-next.best[[minAIC]][1,]
dimnames(top.models)[[3]][z]<-AICS[minAIC]

next.best[[minAIC]]<-next.best[[minAIC]][-1,,drop=F]

if(all(unlist(lapply(next.best,nrow))==0)) {
top.models<-top.models[,,1:z,drop=F]
break}

}}			#  ========= END n top model loop =========


top.models.all<-top.models[,,c(1,
order(as.numeric(dimnames(top.models)[[3]]))[
-which( order(as.numeric(dimnames(top.models)[[3]])) ==1) ]),drop=F]

if(ntop>z) ntop<-z
top.models<-top.models.all[,,1:ntop,drop=F]
class(top.models)<-"MARtop"

top.models

}

