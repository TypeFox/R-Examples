MLEcov<-function(S,nc){
	R=0
	n=sum(nc)
	for(i in 1:length(S)){
	R=R+nc[i]*S[[i]]
	}
	return(R/n)
}

RidgeSample<-function(S,nc,lam){
	p=dim(S[[1]])[1]
	L=eigen(MLEcov(S,nc),symmetric=TRUE)
	NewEVS=(-L$val+sqrt(L$val^2+4*lam))/(2*lam)
	NewOmeg=tcrossprod(L$vec*rep(NewEVS,each=p),L$vec)
	return(NewOmeg)
}

RidgeSampleL<-function(Z,lam){
	S2=lapply(Z,function(x){return(x[[1]])})
	nc2=unlist(lapply(Z,function(x){return(x[[2]])}))
	return(RidgeSample(S2,nc2,lam))
}

RidgeFusion<-setClass("RidgeFusion",slots=c(Omega="list", Ridge="numeric",FusedRidge="numeric",iter="numeric"))

print.RidgeFusion<-function(x,...){
	cat("Number of iterations to convergence",x$iter,"\n")
	cat("Ridge tuning parameter=",x$Ridge,"\n")
	cat("Ridge Fusion tuning parameter=",x$FusedRidge,"\n")
	cat("Precision Matrix Estimates are store in Omega")
}

setMethod("print","RidgeFusion",print.RidgeFusion)

RidgeFused<-function(S,lambda1,lambda2,nc,tol=10^(-7),maxiter=1e3,warm.start=NULL,scale=FALSE){
	iter=0
	J=length(S)
	p=dim(S[[1]])[1]
	n=sum(nc)
	if(is.null(warm.start)){
	Omeg1=lapply(S,function(x){return((diag(1/diag(x),p,p)))})
	}else{Omeg1=warm.start}
	
	
	getsumabs=function(om){return(sum(abs(diag(om))))}
	tolb=tol*sum(sapply(Omeg1,getsumabs))
	
	
	
	if(lambda2>=10^4){
		Init=RidgeSample(S,nc,lambda1)
		for(i in 1:J){
			Omeg1[[i]]=Init
		}
	}
	
		

	Diffiter=1
			Z=0
		for(m in 1:J){
				Z=Z+Omeg1[[m]]
			}
	D=lapply(S,function(x){return(sqrt(1/diag(x)))})
	if(scale==TRUE){S1=lapply(S,cov2cor)}else{S1=S}		
	if(lambda2!=Inf){
	while(Diffiter>tolb && iter<=maxiter ){
	Diffiter=0
	for(j in 1:J){
		Omeg=Omeg1
			lam=(lambda1+lambda2*(J-1))/nc[j]
			L2=eigen(S1[[j]]-(lambda2/nc[j])*(Z-Omeg1[[j]]),symmetric=TRUE)
			NewEVS=(-L2$val+sqrt(L2$val^2+4*lam))/(2*lam)
			NewOmeg=tcrossprod(L2$vec*rep(NewEVS,each=p),L2$vec)
			Diffiter=Diffiter+sum(abs(NewOmeg-Omeg1[[j]]))
			Z=Z-Omeg1[[j]]+NewOmeg
			Omeg1[[j]]=NewOmeg
	##print(Diffiter)
	}
	iter=iter+1
}	

Omeg2=list(0,0)

if(scale==TRUE){
	for(k in 1:J){
		Omeg2[[k]]=rep(D[[k]],each=p)*Omeg1[[k]]*rep(D[[k]],each=p)
	}
}else{Omeg2=Omeg1}
}	
if(lambda2==Inf){
		Omeg1=RidgeSample(S1,nc,lambda1)
		Omeg2=list(0,0)
		if(scale==TRUE){
			
			for(k in 1:J){Omeg2[[k]]=rep(D[[k]],each=p)*Omeg1*rep(D[[k]],each=p)}}else{
				for(k in 1:J){Omeg2[[k]]=Omeg1}}
	}

	##print(iter)
	Rest=list(Omega=Omeg2,Ridge=lambda1,FusedRidge=lambda2,iter=iter)
	class(Rest)="RidgeFusion"	
	return(Rest)
	}
	
 
RidgeFusedL<-function(Z,lambda1,lambda2,tole,ws=NULL,scaleL){
	S2=lapply(Z,function(x){x[[1]]})
	nc2=unlist(lapply(Z,function(x){x[[2]]}))
	Ret=RidgeFused(S2,lambda1,lambda2,nc2,tol=tole,warm.start=ws,scale=scaleL)$Omega
	return(Ret)
}

RidgeFusionCV<-setClass("RidgeFusionCV", slots=c(BestRidge="numeric",BestFusedRidge="numeric",CV="matrix"))


print.RidgeFusionCV<-function(x,...){
	cat("The optimal tuning parameters are Ridge:",x$BestRidge)
	cat("\n Ridge Fusion:",x$BestFusedRidge,"\n")
	cat("The Validation Likelihood Scores\n")
	print(x$CV)
}

setMethod("print","RidgeFusionCV",print.RidgeFusionCV)
	
RidgeFusedCV<-function(X,lambda1,lambda2,Fold,tol=10^-6,warm.start=TRUE,scaleCV=FALSE,INF=FALSE){
	
	S2=lapply(X,function(x){(dim(x)-1)*cov(x)/dim(x)})
	L1=length(lambda1)
	L2=length(lambda2)
	nc=lapply(X,function(x){dim(x)[1]})
	J=length(X)
	nc=lapply(X,function(x){dim(x)[1]})
	p=dim(X[[1]])[2]
	if(length(Fold)==1){stop("Use RidgeFusedRidgeJOICE when only looking for a Single 1.  To use CV Fold>1")}
	
	if(length(Fold)>1){
	K=length(Fold)
	XValid=list(0,0)
	XTest=list(0,0)
	S=XValid
	TV2=XValid	
	
		for(k in 1:K){
			XValid[[k]]=list(0,0)
			XTest[[k]]=list(0,0)
			S[[k]]=list(0,0)
			for(j in 1:J){
				XValid[[k]][[j]]=X[[j]][-Fold[[k]][[j]],]
				XTest[[k]][[j]]=X[[j]][Fold[[k]][[j]],]
			}
			S[[k]]=lapply(XValid[[k]],function(x){return(list((dim(x)[1]-1)*cov(x)/(dim(x)[1]),dim(x)[1]))})
			TV2[[k]]=lapply(XTest[[k]],function(x){(dim(x)[1]-1)*cov(x)/(dim(x)[1])})
		}
		
	Inits=vector("list",L1*L2)
	for(i in 1:L1){
		for(j in 1:L2){
			Inits[[L2*(i-1)+j]]=0
	Inits[[L2*(i-1)+j]]=RidgeFusedL(S[[k]],lambda1[i],Inf,##lambda2[j],
	tole=tol,ws=NULL,scaleL=scaleCV)
	}
	}
		Rest=matrix(0,length(lambda1),length(lambda2))
	for(m in 1:length(lambda1)){
		OmegEst=0
		for(l in 1:length(lambda2)){
			M=0
	OmegEst=list(0,0)		
	OmegEst2=OmegEst
for(k in 1:K){		
	if(warm.start==TRUE){
	OmegEst[[k]]=RidgeFusedL(S[[k]],lambda1[m],lambda2[l],tole=tol,ws=Inits[[L2*(m-1)+l]],
	scaleL=scaleCV)}else{
		OmegEst[[k]]=RidgeFusedL(S[[k]],lambda1[m],lambda2[l],tole=tol,ws=NULL,
	scaleL=scaleCV)
		
		}

				for(j in 1:J){
			M=M+(length(Fold[[k]][[j]])/2)*(sum(diag(TV2[[k]][[j]]%*%OmegEst[[k]][[j]]))-determinant(OmegEst[[k]][[j]])$mod)
				}
			}
		Rest[m,l]=M/K
		}
	}
	if(INF==TRUE){
	InfRest=rep(0,length(lambda1))
	for(m in 1:length(lambda1)){
		M=0
	InfEst=mapply(RidgeSampleL,S,lambda1[m],SIMPLIFY=FALSE)
	for(k in 1:k){
		for(j in 1:J){
			M=M+(length(Fold[[k]][[j]])/2)*(sum(diag(TV2[[k]][[j]]%*%InfEst[[k]]))-log(det(InfEst[[k]])))
		}
	}
	InfRest[m]=M/K
	}
	Rest=cbind(Rest,InfRest)
	lambda2=c(lambda2,Inf)
	}
	colnames(Rest)=lambda2
	row.names(Rest)=lambda1
	Min=which(Rest==min(Rest),arr.ind=TRUE)[1,]
	Lam1=lambda1[Min[1]]
	Lam2=lambda2[Min[2]]
	}
	
	Ret=list(BestRidge=Lam1,BestFusedRidge=Lam2,CV=Rest)
	class(Ret)="RidgeFusionCV"
	return(Ret)
}

RidgeFusedQDA<-setClass("RidgeFusedQDA",slots=c(Omega="list",Means="list",Pi="vector",lambda1="numeric",lambda2="numeric"))	

predict.RidgeFusedQDA<-function(object,newdata,class=TRUE,...){
	x=object
	K=length(x$Omeg)
	n=dim(newdata)[1]
	Result=matrix(nrow=n,ncol=K)
	Class=rep(0,n)
	for(i in 1:n){
		for(k in 1:K){
			Result[i,k]=-.5*(t(newdata[i,]-x$Means[[k]])%*%x$Omeg[[k]]%*%(newdata[i,]-x$Means[[k]])-determinant(x$Omeg[[k]])$mod)+log(x$Pi[[k]])
		}
	Class[i]=which(Result[i,]==max(Result[i,]))	
	}
	
	if(class==TRUE){
		return(list(Predicted=Class))
	}
	if(class==FALSE){
		return(list(Predicted=Result))
	}
}

setMethod("predict","RidgeFusedQDA",predict.RidgeFusedQDA)

print.RidgeFusedQDA<-function(x,...){
cat("The resulting Ridge Fused Quadratic Discriminant Analysis with tuning parameters:\n")
cat("Ridge:",x$lambda1,"\n")
cat("Ridge Fused:",x$lambda2,"\n")
cat("Precision matrices are contained in Omega \n")
cat("Class Proportions\n",x$Pi,"\n")
cat("Means:\n")
print(x$Means)
}
setMethod("print","RidgeFusedQDA",print.RidgeFusedQDA)

FusedQDA<-function(X,Lambda1,Lambda2,scaleC=FALSE){
	
	Means=lapply(X,function(x){apply(x,2,mean)})
	SampCov=lapply(X,function(x){(dim(x)-1)[1]*cov(x)/(dim(x)[1])})
	nc=lapply(X,function(x){dim(x)[1]})
	pi=unlist(nc)/(sum(unlist(nc)))
	Precision=RidgeFused(SampCov,lambda1=Lambda1,lambda2=Lambda2,nc=unlist(nc),scale=scaleC)$Omeg
	Ret=list(Omeg=Precision, Means=Means, Pi=pi,lambda1=Lambda1,lambda2=Lambda2)
	class(Ret)="RidgeFusedQDA"
	return(Ret)	
}
	

logxpy <- function(lp){
	mt=which(lp==max(lp))
	  Ret=lp[mt]
	  for(i in 1:length(lp)){
	  if(i != mt){
	  	Ret=Ret+log1p(exp(-abs(lp[mt]-lp[i])))
	  }
	  }
	  return(Ret)
}

SSRidgeFusion<-setClass("SSRidgeFusion",slots=c(Alphas="matrix",Means="list",Pi="vector"),contains="RidgeFusion")	
print.SSRidgeFusion<-function(x,...){
	cat("The EM Algorithm Converged in ",x$iter,"iterations\n")
	cat("Ridge tuning parameter=",x$Ridge,"\n")
	cat("Ridge Fusion tuning parameter=",x$FusedRidge,"\n")
}

setMethod("print","SSRidgeFusion",print.SSRidgeFusion)

SSRidgeFused<-function(Z,Xu,lambda1,lambda2,Scale=FALSE,warm=NULL,tol=.001){
	J=length(Z)
	p=dim(Z[[1]])[1]
	nc=unlist(lapply(Z,function(x){dim(x)[1]}))
	for(j in 1:J){
		Z[[j]]=as.matrix(Z[[j]])
	}
		Xu=as.matrix(Xu)
	if(is.null(warm)){
	Est=lapply(Z,function(x){(dim(x)[1]-1)*cov(x)/dim(x)[1]})
	Out=RidgeFused(Est,lambda1,lambda2,nc,scale=Scale)$Omega
	Inv=Out
	MeanX1=lapply(Z,function(x){apply(x,2,mean)})
	MeanX=MeanX1
	pi=unlist(nc)/sum(unlist(nc))
	}else{
		Inv=warm$Omega
		MeanX1=warm$Means
		MeanX=MeanX1
		pi=warm$Pi
		}
	alpha=matrix(0,dim(Xu),J)	
    alpha1=alpha+1
    iter=1
    
    ##MeanX=list(0,0)
    Addition=list(0,0)
    Estimate=list(0,0)
    p=dim(Xu)[2]
    
while(max(abs(alpha1-alpha))>=tol && iter<=1e2){

			
			alpha1=alpha	
			Xreas=vector("list",J)

			for(j in 1:J){
			for(v in 1:dim(Xu)[1]){
            Xreas[[j]][v]=(Xu[v,]-MeanX[[j]])%*%Inv[[j]]%*%(Xu[v,]-MeanX[[j]])-determinant(Inv[[j]])$mod
            }
            }
			Tm=list(0,0,0)
		for(v in 1:dim(Xu)[1]){
			Tm[[v]]=rep(0,J)
			for(j in 1:J){
			Tm[[v]][j]=log(pi[j])-.5*Xreas[[j]][v]
			}
		}
		Tm2=lapply(Tm,logxpy)
	for(v in 1:dim(Xu)[1]){
		for(j in 1:J){
			alpha[v,j]=exp(Tm[[v]][j]-Tm2[[v]])
		}
	}

##print("Alpha Calced")

	for(j in 1:J){
		pi[j]=(nc[[j]]+sum(alpha[,j]))/(sum(unlist(nc))+dim(Xu)[1])
	}
   Addition=vector("list",J)
   CX=vector("list",J)
      for(j in 1:J){
   	Addition[[j]]=0
   	CX[[j]]=0
   }
   for(j in 1:J){ 
MeanX[[j]]=(MeanX1[[j]]*nc[[j]]+apply(alpha[,j]*Xu,2,sum))/(nc[[j]]+sum(alpha[,j]))
for(i in 1:dim(Xu)[1]){
Addition[[j]]=Addition[[j]]+alpha[i,j]*(Xu[i,]-MeanX[[j]])%*%t(Xu[i,]-MeanX[[j]])
}
for(i in 1:nc[[j]]){
CX[[j]]=CX[[j]]+(Z[[j]][i,]-MeanX[[j]])%*%t(Z[[j]][i,]-MeanX[[j]])
}
Estimate[[j]]=(Addition[[j]]+CX[[j]])/(nc[[j]]+sum(alpha[,j]))

}

Bottom=list(0,0)
for(j in 1:J){
Bottom[[j]]=nc[[j]]+sum(alpha[,j])
}



Inv=RidgeFused(Estimate,lambda1,lambda2,unlist(Bottom),scale=Scale)$Omega


iter=iter+1
}
##print(iter)

BACK=list(Omega=Inv,Ridge=lambda1,FusedRidge=lambda2, iter=iter, Alpha=alpha,Means=MeanX,Pi=pi)
class(BACK)="SSRidgeFusion"
return(BACK)
}

SSRidgeFusedL<-function(V,lambda1,lambda2,Warm,tolL,scaleL){
	M=SSRidgeFused(V[[1]],V[[2]],lambda1,lambda2,Scale=scaleL,tol=tolL,warm=Warm)
	return(M)
}

predict.SSRidgeFusion<-function(object,newdata,class=TRUE,...){
	x=object
	Look=list(Omega=x$Omega,Means=x$Means,Pi=x$Pi,lambda1=x$Ridge,lambda2=x$FusedRidge)
	class(Look)="RidgeFusedQDA"
	predict(Look,newdata,class)
}

setMethod("predict","SSRidgeFusion",predict.SSRidgeFusion)

SSRidgeFusedCV<-function(X,Xu,Lam1,Lam2,Fold,FoldU,scaleCV=FALSE,tolCV=.01){	
	L1=length(Lam1)
	L2=length(Lam2)
	Inits=vector("list",length(Lam1)*length(Lam2))
	for(i in 1:length(Lam1)){
	for(j in 1:length(Lam2)){
			LH=SSRidgeFused(X,Xu,Lam1[i],Lam2[j],Scale=scaleCV)
			Inits[[L2*(i-1)+j]]=list(Omega=LH$Omega,Ridge=LH$Ridge,FusedRidge=LH$FusedRidge,iter=LH$iter,Alphas=LH$Alpha,Means=LH$Means,Pi=LH$Pi)
		}
	}
	
	Norm=matrix(0,length(Lam1),length(Lam2))
	nc=unlist(lapply(X,function(x){dim(x)[1]}))
	K=length(Fold)
	XValid=list(0,0)
	XTest=list(0,0)
	J=length(X)
	S=XValid
	TV=XValid
	TV2=TV
	
		for(k in 1:K){
			XValid[[k]]=list(0,0)
			XTest[[k]]=list(0,0)
			S[[k]]=list(0,0)
			for(j in 1:J){
				XValid[[k]][[j]]=as.matrix(X[[j]][Fold[[k]][[j]],])
				XTest[[k]][[j]]=as.matrix(X[[j]][Fold[[k]][[j]],])
			}
		}
		
	XTrainU=list(0,0)
	XTestU=list(0,0)
	NU=list(0,0)
	 MV=list(0,0)
     for(k in 1:K){
    MV[[k]]=list(0,0,NULL)
	XTrainU[[k]]=Xu[-FoldU[[k]],]
	XTestU[[k]]=Xu[FoldU[[k]],]
	MV[[k]][[1]]=XValid[[k]]
	MV[[k]][[2]]=XTrainU[[k]]
	}
   
   for(m in 1:length(Lam1)){
   	for(l in 1:length(Lam2)){
   		One=0
   		Two=0
		Here=rep(0,K)
		Mod=lapply(MV,SSRidgeFusedL,Lam1[m],Lam2[l],scaleL=scaleCV,tolL=tolCV,Warm=Inits[[L2*(m-1)+l]])
		## Put Inits Here 

   		for(k in 1:K){
   			MV[[k]][[3]]=Mod[[k]]$Alpha
   			TM=list(0,0)
   			AlphU=matrix(0,dim(XTestU[[k]])[1],J)	
   			for(i in 1:dim(XTestU[[k]])[1]){
   				for(j in 1:J){
   				AlphU[i,j]=-log(Mod[[k]]$Pi[[j]])+.5*((XTestU[[k]][i,]-Mod[[k]]$Mean[[j]])%*%Mod[[k]]$Omega[[j]]%*%(XTestU[[k]][i,]-Mod[[k]]$Mean[[j]])-determinant(Mod[[k]]$Omega[[j]])$mod[1])
   			}
   			##One=One+log(sum(AlphU[i,]))
   		}
   		One=One+sum(apply(AlphU,1,logxpy))
   		Here[k]=sum(apply(AlphU,1,logxpy))
   		for(j in 1:J){
   			for(d in 1:dim(XTest[[k]][[j]])[1]){
   				
   				Two=Two-log(Mod[[k]]$Pi[[j]])+.5*((XTest[[k]][[j]][d,]-Mod[[k]]$Means[[j]])%*%Mod[[k]]$Omega[[j]]%*%(XTest[[k]][[j]][d,]-Mod[[k]]$Means[[j]])-determinant(Mod[[k]]$Omega[[j]])$mod[1])
   			   			}
   		}
   		}
   		Norm[m,l]=(1/K)*(-One+Two)
   	}
   }
Min=which(Norm==min(Norm),arr.ind=TRUE)
	colnames(Norm)=Lam2
	rownames(Norm)=Lam1
	L1=Lam1[Min[1]]
	L2=Lam2[Min[2]]
	
	Ret=list(BestRidge=L1,BestFusedRidge=L2,CV=Norm)
    class(Ret)="RidgeFusionCV"
    return(Ret)
}