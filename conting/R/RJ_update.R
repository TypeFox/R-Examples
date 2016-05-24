RJ_update <-
function(prop.index,curr.index,curr.beta,eta.hat,curr.y,big.X,proposal.probs,i.prop.prior.var,i.curr.prior.var){

icurrR<-i.curr.prior.var[-1,-1] 				## current inverse prior variance
ipropR<-i.prop.prior.var[-1,-1]					## proposed inverse prior variance
currR<-chol2inv(chol(icurrR))					## current prior variance
propR<-chol2inv(chol(ipropR))					## proposed prior variance

RHO_TOP<-proposal.probs[1]					## proposal prob (proposed to current)
RHO_BOT<-proposal.probs[2]					## proposal prob (current to proposed)

curr.X<-big.X[,curr.index==1]					## current design matrix
prop.X<-big.X[,prop.index==1]					## proposed design matrix
S<-matrix(big.X[,prop.index!=curr.index],nrow=length(curr.y))	## difference in current and proposed design matrices
w<-exp(eta.hat)							## weight matrix under maximal model

curr.LP<-as.vector(curr.X%*%matrix(curr.beta,ncol=1))		## current linear predictor

if(sum(prop.index)<sum(curr.index)){  ### Death move

pXW<-t(prop.X)%*%diag(w)
pXWX<-pXW%*%prop.X
ipXWX<-chol2inv(chol(pXWX))
ipXWX.XW<-ipXWX%*%pXW
SWIP<-t(S)%*%diag(w)-t(S)%*%t(pXW)%*%ipXWX.XW
SIG1<-chol2inv(chol(SWIP%*%S))					## proposal variance
MU1<-as.vector(SIG1%*%SWIP%*%matrix(eta.hat,ncol=1))		## proposal mean

beta_1<-curr.beta[prop.index[curr.index==1]==1]			## parameters to retain
beta_2<-curr.beta[prop.index[curr.index==1]==0]			## parameters to kill

prop.beta<-beta_1+as.vector(ipXWX.XW%*%S%*%matrix(beta_2,ncol=1))	## proposal
prop.LP<-as.vector(prop.X%*%matrix(prop.beta,ncol=1))			## proposed linear predictor

top<-sum(curr.y*prop.LP)-sum(exp(prop.LP))+dmvnorm(x=prop.beta[-1],mean=rep(0,length(prop.beta)-1),sigma=propR,log=TRUE)
bot<-sum(curr.y*curr.LP)-sum(exp(curr.LP))+dmvnorm(x=curr.beta[-1],mean=rep(0,length(curr.beta)-1),sigma=currR,log=TRUE)
jac<-dmvnorm(x=beta_2,mean=MU1,sigma=SIG1,log=TRUE) 		## log numerator, denominator and jacobian

prob<-(RHO_TOP/RHO_BOT)*exp(top-bot+jac)}			## acceptance probability

if(sum(prop.index)>sum(curr.index)){  ### Birth move

cXW<-t(curr.X)%*%diag(w)
cXWX<-cXW%*%curr.X
icXWX<-chol2inv(chol(cXWX))
icXWX.XW<-icXWX%*%cXW
SWIP<-t(S)%*%diag(w)-t(S)%*%t(cXW)%*%icXWX.XW
SIG1<-chol2inv(chol(SWIP%*%S))					## proposal variance
MU1<-as.vector(SIG1%*%SWIP%*%matrix(eta.hat,ncol=1))		## proposal mean	

u1<-as.vector(rmvnorm(n=1,mean=MU1,sigma=SIG1))			## innovation variables

prop.beta<-rep(0,dim(prop.X)[2])
prop.beta[curr.index[prop.index==1]==1]<-curr.beta-as.vector(icXWX.XW%*%S%*%matrix(u1,ncol=1))
prop.beta[curr.index[prop.index==1]==0]<-u1			## proposal for beta
prop.LP<-as.vector(prop.X%*%matrix(prop.beta,ncol=1))		## proposed linear predictor

top<-sum(curr.y*prop.LP)-sum(exp(prop.LP))+dmvnorm(x=prop.beta[-1],mean=rep(0,length(prop.beta)-1),sigma=propR,log=TRUE)
bot<-sum(curr.y*curr.LP)-sum(exp(curr.LP))+dmvnorm(x=curr.beta[-1],mean=rep(0,length(curr.beta)-1),sigma=currR,log=TRUE)
jac<--dmvnorm(x=u1,mean=MU1,sigma=SIG1,log=TRUE)		## log numerator, denominator and jacobian

prob<-(RHO_TOP/RHO_BOT)*exp(top-bot+jac)}			## acceptance probability

if(prob>=runif(1)){
new.beta<-prop.beta
new.index<-prop.index} else{
new.beta<-curr.beta
new.index<-curr.index}						## accept or reject

list(new.beta=new.beta,new.index=new.index)}
