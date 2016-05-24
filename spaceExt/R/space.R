
######implementation of glasso with EM
glasso.miss<-function(Y, emIter, rho, ...){
#zero=NULL, thr=1.0e-4, maxit=1e4,  approx=FALSE, penalize.diagonal=TRUE, start=c("cold","warm"), w.init=NULL,wi.init=NULL, trace=FALSE
##zero, thr, maxit,  approx, penalize.diagonal, start, w.init,wi.init, trace
missInfo=is.na(Y)
n=nrow(Y)
p=ncol(Y)
if(sum(missInfo)==0) emIter=0
Y[missInfo]=0
s=matrix(0,p,p)
lastOmega=matrix(0,p,p)
i=0
while(i < emIter){
i=i+1
print(paste("Em step:",i))
s=s/(n)+(n-1)/n*cov(Y)
result=glasso(s,rho, ...)
s=matrix(0,p,p)	
	Omega=result$wi
	for(emE in 1:n){ #expecation step
		emMiss<-missInfo[emE,]
		if(sum(emMiss)>1){
			emOmega<-Omega[emMiss,emMiss]
			condiVar<-Solve(emOmega)
			s[emMiss,emMiss]<-condiVar+s[emMiss,emMiss]
			emB<-Omega[emMiss,!emMiss]
			Y[emE,emMiss]=-(condiVar%*%emB%*%Y[emE,!emMiss])
		}else if(sum(emMiss)>0){
			emOmega<-Omega[emMiss,emMiss]
			emB<-Omega[emMiss,!emMiss]
			s[emMiss,emMiss]<-1/emOmega+s[emMiss,emMiss]
			Y[emE,emMiss]=(-emB%*%Y[emE,!emMiss])/emOmega
		}
	} ###end of expectation loop
maxDif<-max(abs(Omega-lastOmega))
	if(maxDif<0.001){
		print("Em Converge!")
		i=emIter
	}
lastOmega=Omega
}###end of em iteration


##bic

	temp.s=0
  	for(emE in 1:n){ #expecation step
		emMiss<-missInfo[emE,]
		Gamma=Omega[!emMiss,!emMiss]
		if(sum(emMiss)>1){
			emOmega<-Omega[emMiss,emMiss]
			condiVar<-Solve(emOmega)
			Gamma=Gamma-Omega[!emMiss,emMiss]%*%condiVar%*%Omega[emMiss,!emMiss]
		}else if(sum(emMiss)>0){
			emOmega<-Omega[emMiss,emMiss]
			emB<-as.matrix(Omega[emMiss,!emMiss],ncol=1)
			Gamma=Gamma-emB%*%t(emB)/emOmega
		}
		temp.s=temp.s-log(det(Gamma))+Y[emE,!emMiss]%*%Gamma%*%Y[emE,!emMiss]
	}
  Omega[(abs(Omega)<0.001)]=0;diag(Omega)=0
  temp.s<-temp.s+log(n)*(sum(Omega!=0)/2+p)
result=c(list(Y.imputed=Y,bic=temp.s),result)
return(result)
}




space.miss<-function(Y.m, lam1, lam2=0, sig=NULL, weight=NULL,iter=2,emIter=5,n_iter=1000,t0=0,r=0)
{
  n=nrow(Y.m)
  p=ncol(Y.m)
  
  
ITER=iter
################### preparation###########
if(!is.null(sig))
  { #### if the user specify "sig", do not need to update sig.
     SIG.update=F
     SIG=sig
  } else { #### otherwise, need to update sig in each iteration
     SIG.update=T
     SIG=rep(1, p) 
  } 
  
if(length(weight)==0 | (length(weight)>1 & length(weight)<p)) ### weight==NULL 
  {
    WEIGHT.tag=0 ### equal weight 
    WEIGHT=rep(1, p)
    WEIGHT.update=F
  } 
if(length(weight)==1)
  {
    if(weight==1)
     {
      WEIGHT.tag=1 ### sig based weight
      WEIGHT=SIG
      WEIGHT.update=T
      ITER=max(2, iter)
     } else {
      WEIGHT.tag=2 ### degree based weight
      WEIGHT=rep(1,p)
      WEIGHT.update=T
      ITER=max(2, iter)
     } 
  }
if(length(weight)==p)  
  {
     WEIGHT.tag=3 ## prespecified weight
     WEIGHT=weight
     WEIGHT.update=F
  }
  
  missInfo<-is.na(Y.m)
  if(!any(missInfo)){
  emIter=0
  }else
  Y.m[missInfo]=0
  
  smoothWeight<-exp(-r*abs(1:n-t0))
  Y.m<-diag(sqrt(smoothWeight))%*%Y.m
  
####initial estimate######
################## begin to iterate
for(i in 1:ITER)
  {    
    Y.u<-Y.m*matrix(sqrt(WEIGHT),n,p,byrow=TRUE)
    sig.u<-SIG/WEIGHT
	
    ParCor.fit<-jsrm(Y.u,missInfo,sig.u,n,p,lam1,lam2,n_iter=n_iter,haction=NULL,Beta_output=NULL)
    diag(ParCor.fit)<-1
    
    coef<-ParCor.fit[upper.tri(ParCor.fit)]
    beta.cur<-Beta.coef(coef,SIG) 
 
  
    if(!WEIGHT.update & !SIG.update)
     {
       break
     } else { ## update sig or weight
        if(SIG.update)
         {
             SIG<-InvSig.diag.new(Y.m,beta.cur,missInfo)   
         }
        if(WEIGHT.update)
         {
             if(WEIGHT.tag==1)
             {        #### sig based
               WEIGHT=SIG
             } else { #### degree based
               temp.w<-apply(abs(ParCor.fit)>1e-6,1,sum)
               temp.w<-temp.w+max(temp.w)     
               WEIGHT<-temp.w/sum(temp.w)*p    
             }
          } 
     } ### end updating WEIGHT and SIG
  } ### end iteration

  result<-list(ParCor=ParCor.fit,sig.fit=SIG,weight.fit=WEIGHT)
  lastParCor=ParCor.fit
  
#####beginning of EM algorithm
  em=1
	while(em <=emIter){
	print(paste("EM steps:",em))
	hterm=matrix(0,p,p)  ##additional term
	sig.fit<-result$sig.fit;weight.fit<-result$weight.fit
	temp=-(result$ParCor);diag(temp)=1;temp[(abs(temp)<0.0001)]=0
	Omega=diag(sqrt(sig.fit))%*%temp%*%diag(sqrt(sig.fit))
	
	for(emE in 1:n){ #expecation step
		emMiss<-missInfo[emE,]
		if(sum(emMiss)>1){
			emOmega<-Omega[emMiss,emMiss]
			condiVar<-Solve(emOmega)
			hterm[emMiss,emMiss]=hterm[emMiss,emMiss]+smoothWeight[emE]*condiVar
			emB<-Omega[emMiss,!emMiss]
			Y.m[emE,emMiss]=-(condiVar%*%emB%*%Y.m[emE,!emMiss])
		}else if(sum(emMiss)>0){
			emOmega<-Omega[emMiss,emMiss]
			emB<-Omega[emMiss,!emMiss]
			hterm[emMiss,emMiss]=hterm[emMiss,emMiss]+smoothWeight[emE]/emOmega
			Y.m[emE,emMiss]=(-emB%*%Y.m[emE,!emMiss])/emOmega
		}
	}
	#Maximization Step
	hterm=diag(sqrt(sig.fit))%*%hterm%*%diag(sqrt(sig.fit))
	result=space.max(Y.m,lam1=lam1, lam2=lam2,haction=hterm,n_iter=n_iter, parCor=lastParCor,
					sig=sig.fit,SIG.update=SIG.update,weight=weight.fit,WEIGHT.tag=WEIGHT.tag,WEIGHT.update=WEIGHT.update) 
	maxDif<-max(abs(result$ParCor-lastParCor))
	if(maxDif<0.001){
		print("Em Converge!")
		em=emIter
	}
	lastParCor=result$ParCor

	em=em+1
	}##end of EM

	
####BIC####
	sig.fit<-result$sig.fit;

	temp=-(result$ParCor);diag(temp)=1;temp[(abs(temp)<0.0001)]=0
	Omega=diag(sqrt(sig.fit))%*%temp%*%diag(sqrt(sig.fit))
	temp.s=0
  	for(emE in 1:n){ #expecation step
		emMiss<-missInfo[emE,]
		Gamma=Omega[!emMiss,!emMiss]
		if(sum(emMiss)>1){
			emOmega<-Omega[emMiss,emMiss]
			condiVar<-Solve(emOmega)
			Gamma=Gamma-Omega[!emMiss,emMiss]%*%condiVar%*%Omega[emMiss,!emMiss]
		}else if(sum(emMiss)>0){
			emOmega<-Omega[emMiss,emMiss]
			emB<-as.matrix(Omega[emMiss,!emMiss],ncol=1)
			Gamma=Gamma-emB%*%t(emB)/emOmega
		}
		temp.s=temp.s-smoothWeight[emE]*log(det(Gamma))+Y.m[emE,!emMiss]%*%Gamma%*%Y.m[emE,!emMiss]
	}
  temp[(abs(temp)<0.001)]=0;diag(temp)=0
  temp.s<-temp.s+log(n)*(sum(temp!=0)/2+p)
  result<-list(Y.imputed=Y.m,ParCor=result$ParCor,sig.fit=result$sig.fit,bic=temp.s,
				SIG.update=SIG.update,weight=weight.fit,WEIGHT.tag=WEIGHT.tag,WEIGHT.update=WEIGHT.update,
				n_iter=n_iter)
  return(result)  
}


space.max<-function(Y.m, lam1, lam2=0,haction=NULL,n_iter=1000,parCor=NULL,sig,SIG.update,weight,WEIGHT.tag,WEIGHT.update)
{
SIG=sig
WEIGHT=weight

n=nrow(Y.m)
p=ncol(Y.m)

    Y.u<-Y.m*matrix(sqrt(WEIGHT),n,p,byrow=TRUE)
    sig.u<-SIG/WEIGHT
	
    ParCor.fit<-jsrm(Y.u,missInfo=NULL,sig.u,n,p,lam1,lam2,n_iter=n_iter,haction=haction,Beta_output=parCor)
    diag(ParCor.fit)<-1
    
    coef<-ParCor.fit[upper.tri(ParCor.fit)]
    beta.cur<-Beta.coef(coef,SIG) 
 
  
    if(!WEIGHT.update & !SIG.update)
     {
       break
     } else { ## update sig or weight
        if(SIG.update)
         {
             SIG<-InvSig.diag.new(Y.m,beta.cur,missInfo=NULL)   
         }
        if(WEIGHT.update)
         {
             if(WEIGHT.tag==1)
             {        #### sig based
               WEIGHT=SIG
             } else { #### degree based
               temp.w<-apply(abs(ParCor.fit)>1e-6,1,sum)
               temp.w<-temp.w+max(temp.w)     
               WEIGHT<-temp.w/sum(temp.w)*p    
             }
          } 
     } ### end updating WEIGHT and SIG

  result<-list(ParCor=ParCor.fit,sig.fit=SIG,weight.fit=WEIGHT)
  return(result)  
}



########################
##################################################################
############################ internal functions 
#################################################################

jsrm<-function(Y,missInfo=NULL,sig.use,n,p,lam1,lam2=0, n_iter=1000,haction=NULL,Beta_output=NULL)
{
lambda1=lam1
lambda2=lam2
sigma_sr=sig.use^0.5


if(is.null(Beta_output))
	Beta_output<-rep(0, p*p)
else{
	diag(Beta_output)=0
   Beta_output<-as.vector(t(Beta_output))
   }
   
Y_data<-as.vector(t(Y))

if(!is.null(missInfo))
missInfo_data<-as.vector(t(missInfo))

if(!is.null(haction)) haction_data<-as.vector(t(haction))

#dyn.load("JSRM.so") ### compiled from "JSRM.c"
                   
testConver=0
if(!is.null(haction)){
junk<-.C("JSRM2", 
           as.integer(n),
           as.integer(p),
           as.single(lambda1),
           as.single(lambda2),
           as.single(Y_data),
           as.single(sigma_sr),
           as.integer(n_iter),
           as.single(haction_data),
		   beta.estimate=as.single(Beta_output)
		 )
}else{
junk<-.C("JSRM", 
		   as.integer(missInfo_data),
           as.integer(n),
           as.integer(p),
           as.single(lambda1),
           as.single(lambda2),
           as.single(Y_data),
           as.single(sigma_sr),
           as.integer(n_iter),
           beta.estimate=as.single(Beta_output),
		   testConver=as.single(testConver)
         )
}

#print(testConver)			 
beta.v<-junk$beta.estimate
beta.m<-matrix(beta.v, p,p, byrow=T)
return(beta.m)
}

########################################
##### Estimate diagonal sigma
##use the fact that 1/sig^{ii} is the residual variance in one versus all others setting


InvSig.diag.new<-function(Y, Beta,missInfo){
################### parameters
### Y:    n by p data matrix; should be standardized to sd 1;
### Beta: beta matrix for the regression model
 p=ncol(Y)
 Beta.temp=Beta
 diag(Beta.temp)=0
 esti.Y=Y%*%Beta.temp
 residue=Y-esti.Y
 residue[missInfo]=NA
 result=apply(residue^2,2,mean,na.rm=TRUE)
 return(1/result)
}




########################################
### given rho^{ij}, get beta[j,i]: each column for each variable  
### beta[j,i]<-rho^{ij}sqrt(sig^{jj}/sig{ii})

Beta.coef<-function(coef, sig.fit){
############## parameter
### coef: rho^{ij}; 
### sig.fit: sig^{ii}
 p<-length(sig.fit)
 result<-matrix(0,p,p)
 result[upper.tri(result)]<-coef
 result<-result+t(result)
 result<-diag(1/sqrt(sig.fit))%*%result%*%diag(sqrt(sig.fit))
 result<-t(result) 
 return(result)
}



##################################################
## get covariance matrix from a partial correlation matrix 
## make it p.d.

GenCov<-function(ParCor.m){
temp<-(-ParCor.m)
diag(temp)<-1
Sig<-solve(ParCor.m)


p<-nrow(Sig)
for (i in 2:p){
 for (j in 1:(i-1)){
  Sig[i,j]<-Sig[i,j]/sqrt(Sig[i,i]*Sig[j,j]) ##standardize to sigma=1
  Sig[j,i]<-Sig[i,j]
 }
}
diag(Sig)<-1

#diagonose
D<-eigen(Sig)$values
if(!all(D>0)){
print("invalid covariance matrix")
return()
}

return(Sig)
}

