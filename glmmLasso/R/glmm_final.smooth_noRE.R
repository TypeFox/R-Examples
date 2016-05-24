glmm_final_smooth_noRE<-function(y,X,Phi,penal.vec,K,Delta_start,steps=1000,family,overdispersion,
                            phi,nue=1,print.iter.final=FALSE,eps.final=1e-5)
{
dim.smooth<-dim(Phi)[2]  
N<-length(y)
lin<-ncol(as.matrix(X))
Eta<-cbind(X,Phi)%*%Delta_start

if(is.null(family$multivariate)){
  D<-family$mu.eta(Eta)
  Mu<-family$linkinv(Eta)
  SigmaInv <- 1/family$variance(Mu)
}else{
  Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
  Mu_cat <- family$linkinv(Eta_cat)
  D <- family$deriv.mat(Mu_cat)
  SigmaInv <- family$SigmaInv(Mu_cat)
  Mu <- c(t(Mu_cat))
}

if(print.iter.final)
  message("Final Re-estimation Iteration ", 1)

Z_alles<-cbind(X,Phi)

P1<-c(rep(0,lin),penal.vec)
P1<-diag(P1)

Delta<-matrix(0,steps,(lin+dim.smooth))
Eta.ma<-matrix(0,steps+1,N)
Eta.ma[1,]<-Eta

l=1
opt<-steps


if(is.null(family$multivariate)){
  D <- drop(D);SigmaInv <- drop(SigmaInv)
  score_vec <- t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[1,]
  F_gross<-t(Z_alles)%*%(Z_alles*D*SigmaInv*D)+P1
}else{
  score_vec <- t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[1,]
  F_gross<-t(Z_alles)%*%(D%*%(SigmaInv%*%(t(D)%*%Z_alles)))+P1
}

InvFisher<-try(chol2inv(chol(F_gross)),silent=TRUE)
if(class(InvFisher)=="try-error")
InvFisher<-try(solve(F_gross),silent=TRUE)  
if(class(InvFisher)=="try-error")
stop("Fisher matrix not invertible")

half.index<-0
solve.test<-FALSE
Delta_r<-InvFisher%*%score_vec

######### big while loop for testing if the update leads to Fisher matrix which can be inverted
while(!solve.test)
{  
  if(half.index>50)
    half.index<-Inf
  
Delta[1,]<-Delta_start+nue*(0.5^half.index)*Delta_r

Eta<-Z_alles%*%Delta[1,]

if(is.null(family$multivariate)){
  D<-family$mu.eta(Eta)
  Mu<-family$linkinv(Eta)
  SigmaInv <- 1/family$variance(Mu)
}else{
  Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
  Mu_cat <- family$linkinv(Eta_cat)
  D <- family$deriv.mat(Mu_cat)
  SigmaInv <- family$SigmaInv(Mu_cat)
  Mu <- c(t(Mu_cat))
}

if(is.null(family$multivariate)){
  D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
  score_vec <- t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[1,]
  F_gross<-t(Z_alles)%*%(Z_alles*D*SigmaInv*D)+P1
}else{
  score_vec <- t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[1,]
  F_gross<-t(Z_alles)%*%(D%*%(SigmaInv%*%(t(D)%*%Z_alles)))+P1
}

InvFisher<-try(chol2inv(chol(F_gross)),silent=TRUE)
if(class(InvFisher)=="try-error")
  InvFisher<-try(solve(F_gross),silent=TRUE)  
if(class(InvFisher)=="try-error")
{
  half.index<-half.index+1  
}else{
  solve.test<-TRUE 
}
}

Eta.ma[2,]<-Eta

y_dach<-as.vector(family$linkinv(Eta))
Dev_neu<-sum(family$dev.resids(y,y_dach,wt=rep(1,N))^2)

###############################################################################################################################################
################################################################### Main Iteration ###################################################################
eps<-eps.final*sqrt(length(Delta_r))

for (l in 2:steps)
{
if(print.iter.final)
  message("Final Re-estimation Iteration ", l)
#print(paste("Final Re-estimation Iteration ", l,sep=""))

half.index<-0
solve.test<-FALSE


Delta_r<-InvFisher%*%score_vec
######### big while loop for testing if the update leads to Fisher matrix which can be inverted
first.time<-FALSE
while(!solve.test)
{  
  
  if(half.index>50)
    half.index<-Inf
  
Delta[l,]<-Delta[l-1,]+nue*(0.5^half.index)*Delta_r

Eta<-Z_alles%*%Delta[l,]

if(is.null(family$multivariate)){
  D<-family$mu.eta(Eta)
  Mu<-family$linkinv(Eta)
  SigmaInv <- 1/family$variance(Mu)
}else{
  Eta_cat <- matrix(Eta, byrow = TRUE, ncol = K)
  Mu_cat <- family$linkinv(Eta_cat)
  D <- family$deriv.mat(Mu_cat)
  SigmaInv <- family$SigmaInv(Mu_cat)
  Mu <- c(t(Mu_cat))
}

if(is.null(family$multivariate)){
  D <- as.vector(D);SigmaInv <- as.vector(SigmaInv)
  score_vec <- t(Z_alles)%*%((y-Mu)*D*SigmaInv)-P1%*%Delta[l,]
  F_gross<-t(Z_alles)%*%(Z_alles*D*SigmaInv*D)+P1
}else{
  score_vec <- t(Z_alles)%*%(D%*%(SigmaInv%*%(y-Mu)))-P1%*%Delta[l,]
  F_gross<-t(Z_alles)%*%(D%*%(SigmaInv%*%(t(D)%*%Z_alles)))+P1
}

InvFisher<-try(chol2inv(chol(F_gross)),silent=TRUE)
if(class(InvFisher)=="try-error")
  InvFisher<-try(solve(F_gross),silent=TRUE)  
if(class(InvFisher)=="try-error")
{
  half.index<-half.index+1  
}else{
  solve.test<-TRUE 
}
}

Eta.ma[l+1,]<-Eta

finish<-(sqrt(sum((Eta.ma[l,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l,])^2))<eps)
finish2<-(sqrt(sum((Eta.ma[l-1,]-Eta.ma[l+1,])^2))/sqrt(sum((Eta.ma[l-1,])^2))<eps)

if(finish ||  finish2) 
  break
}

######## Final calculation

opt<-l

if(is.null(family$multivariate)){
  W_opt <- D*SigmaInv*D
  FinalHat<-(Z_alles*sqrt(W_opt))%*%InvFisher%*%t(Z_alles*sqrt(W_opt))
}else{
  W_opt <- D%*%(SigmaInv%*%t(D))
  W_inv_t <- chol(W_opt)
  FinalHat<-W_inv_t%*%(Z_alles%*%(InvFisher%*%(t(Z_alles)%*%t(W_inv_t))))
}

complexity<-sum(diag(FinalHat))

if(overdispersion)
  phi<-(sum((y-Mu)^2/family$variance(Mu)))/(N-complexity)


Deltafinal<-Delta[l,]
Standard_errors<-InvFisher

#### compare Wood's Book on page 167
  if(is.null(family$multivariate)){
    EDF.matrix<-InvFisher%*%(t(Z_alles)%*%(Z_alles*D*SigmaInv*D))
  }else{
    EDF.matrix<-InvFisher%*%(t(Z_alles)%*%(D%*%(SigmaInv%*%(t(D)%*%Z_alles))))
  }
ret.obj=list()
ret.obj$opt<-opt
ret.obj$Delta<-Deltafinal
ret.obj$Standard_errors<-Standard_errors
ret.obj$phi<-phi
ret.obj$EDF.matrix<-EDF.matrix
ret.obj$complexity<-complexity
return(ret.obj)
}


