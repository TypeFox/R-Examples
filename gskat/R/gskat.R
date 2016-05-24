# TODO: GSKAT main function block
# 
# Author: Xuefeng Wang
###############################################################################


#parameter estiamtion under the null
gskat_fit<-function(y,XC,ID,F1=FALSE){
	#y: binary phenotype coded as 0, 1
	#XC:covaraite matrix, including PC and intercept 
	#ID data.frame: including four columns FID, IID, PAT , MAT
	#F1: If ture, using identity working corr. matrix
	#OUTPUT: alpha: estimates (under the null); delta; iteraction no.
	#method: Fisher scoring and mom iteration
	
	names(ID)=c("FID","IID","PAT","MAT")
	
	nx<-ncol(XC)#no. of covaraites (including intercpet)	
	#alpha<-t(t(c(-1,rep(1,(ncol(XC)-1))))) #inital values
	alpha<-summary(glm(y~XC[,-1],family=binomial("logit")))$coeff[,1]  #use glm to get initial values	
	d <-0.5
	diff1=diff2=10
	i<-0
	
	while (diff1>10e-6 & diff2>10e-6){		
		i<-i+1		
		#Fisher scoring
		mu<- plogis( XC %*% alpha)
		mu=as.vector(mu)
		mydata<-data.frame(ID, X=XC, y=y, mu=mu)
		result<-Reduce("+",lapply(split(mydata,mydata$FID),function (adata) {
							# adata<-mydata[1:12,]
							yi<-adata[,(ncol(adata)-1)]; mui<-adata[,ncol(adata)]; Xi<-as.matrix(adata[,5:(5+nx-1)])																			
							#Di<-diag((mui)*(1-mui))%*%Xi  #Delta%*%X
							if(length(mui)>1){Di<-diag((mui)*(1-mui))%*%Xi; Ai<-diag((mui)*(1-mui))} else {
								Di<-(mui)*(1-mui)*Xi; Ai<-(mui)*(1-mui)
							}							
							Rd<-mykin(adata$IID,adata$PAT,adata$MAT,F1=F1) ### <-----kinship working correlation
							if(length(Rd)>1){Rd<-Rd*d;diag(Rd)=1}
							Vi<-sqrt(Ai)%*%(Rd)%*%sqrt(Ai)  ###<-------
							Ui<-t(Di)%*%solve(Vi)%*%t(t((yi-mui))) 
							Ii<-t(Di)%*%solve(Vi)%*%Di
							return(cbind(Ui,Ii))
						}))
		U<-t(t(result[,1])); I<-result[,2:ncol(result)]
		rm(mydata,result)
		alpha.new<-alpha+solve(I)%*%U
		diff1 <- max(abs(alpha.new-alpha))
		alpha<-alpha.new
		
		# Method of moments
		mu.new<-plogis(XC%*%alpha)
		r.hat<-(y-mu.new)/sqrt((mu.new*(1-mu.new)))    #Pearson residual
		
		mydata<-data.frame(ID,r.hat=r.hat)		
		result.mom<-Reduce("+",lapply(split(mydata,mydata$FID),function (adata) {
							ri<-adata[,5]
							RI<-tcrossprod(ri)
							mynu<-sum(RI[upper.tri(RI)])
							R1<-mykin(adata$IID,adata$PAT,adata$MAT,F1=F1)
							myde<-sum(R1[upper.tri(R1)])
							return(c(mynu,myde))
						}))
		d.new<-result.mom[1]/result.mom[2]
		diff2 <- abs(d.new - d)		
		d<-d.new
	}	
	return(list(alpha,d,i))
}


#GEE_SKAT score test function
gskat_score<-function(pedDat,F1=FALSE) {
	#pedDat is a list including four data matrix: ID, y, X, Z (see the example dataset)
	#y: binary phenotype coded as 0, 1
	#X:covaraite matrix, including PC and intercept 
	#ID data.frame: including four columns FID, IID, PAT , MAT
	#F1: If ture, using identity working corr. matrix
	ID=pedDat$ID
	y=pedDat$y
	XC=pedDat$X
	Z=pedDat$Z
	
	#require("CompQuadForm")	
	#
	y<-as.vector(y)
	names(ID)=c("FID","IID","PAT","MAT")
	
	#
	if (sum(XC[,1]!=1)!=0) {XC=cbind(1,data.matrix(XC))} #now with intercept, assume all nonvariants removed
	Z<-data.matrix(Z)	#Z is the combined genotype, must be matrix format
	y<-as.vector(y);FID<-as.factor(ID$FID)
	
	
	##sort according to FID
	FID.old<-FID
	FID<-FID[order(FID.old)]; y<-y[order(FID.old)];XC<-XC[order(FID.old),];Z<-data.matrix(Z[order(FID.old),]) #!!
	ID<-ID[order(FID.old),]
	p=ncol(Z); q=ncol(XC)-1; nx=q+1 
	#
	null<-gskat_fit(y,XC,ID,F1=F1)
	
	
	alpha<-null[[1]]; d<-null[[2]]
	mu<-as.vector(plogis(XC%*%alpha))      #mu<-as.vector(plogis(X%*%alpha))
	
	mydata<-data.frame(ID, X=XC, Z=Z, y=y, mu=mu)	
	result<-lapply(split(mydata,mydata$FID),function (adata) {
				#adata=mydata[5,]
				yi<-adata[,(ncol(adata)-1)]; mui<-adata[,ncol(adata)]; Xi<-as.matrix(adata[,5:(5+nx-1)])	
				Zi<-as.matrix(adata[,(ncol(adata)-p-1):(ncol(adata)-2)])
				XZ<-cbind(Xi,Zi)
				
				if(length(mui)>1){Di<-diag((mui)*(1-mui))%*%XZ; Ai<-diag((mui)*(1-mui))} else {
					Di<-(mui)*(1-mui)*XZ; Ai<-(mui)*(1-mui)
				}	
				
				Rd<-mykin(adata$IID,adata$PAT,adata$MAT,F1=F1)
				if(length(Rd)>1){Rd<-Rd*d;diag(Rd)=1}
				Vi<-sqrt(Ai)%*%(Rd)%*%sqrt(Ai)
				iVi<-solve(Vi)
				
				Ui<-t(Zi)%*%Ai%*%iVi%*%t(t((yi-mui)))
				
				Covyi<-tcrossprod(yi-mui)
				#Delta%*%(X,Z)
				Bi<-t(Di)%*%iVi%*%Covyi%*%iVi%*%Di
				AAi<-t(Di)%*%iVi%*%Di
				return(list(Ui,Bi,AAi))
			})
	U=result[[1]][[1]];B=result[[1]][[2]];A=result[[1]][[3]]
	n=length(unique(ID$FID))
	for (i in 2:n) {
		U=U+result[[i]][[1]]
		B=B+result[[i]][[2]]
		A=A+result[[i]][[3]]
	}
	
	Azx<-A[(q+1+1):ncol(A),1:(q+1)]
	Axx<-A[1:(q+1),1:(q+1)]
	C<-cbind(-Azx%*%solve(Axx),diag(p))
	
	TS<-t(U)%*%U
	Bsqrt=mysqrt(B)
	BC<-Bsqrt%*%t(C)%*%C%*%Bsqrt
	
	Lamda<-eigen(BC,only.values=T)$values
	fout<-davies(TS,Lamda,rep(1, length(Lamda)))
	return(list(pval=fout$Qq,ifault=fout$ifault))
}


#########################################
# perturbation test
##########################################
gskat_score_pert<-function(pedDat,F1=FALSE,pw="Rade",np=10000) {
	#pedDat is a list including four data matrix: ID, y, X, Z (see the example dataset)
	#y: binary phenotype coded as 0, 1
	#X:covaraite matrix, including PC and intercept 
	#ID data.frame: including four columns FID, IID, PAT , MAT
	#F1: If ture, using identity working corr. matrix
	#pw: peturbation method: "Rade"=Rademacher, "Norm"=Normal distribution
	#np: number of perturbed samples
	
	#require(e1071)
	#require(Matrix)
	
	ID=pedDat$ID
	y=pedDat$y
	XC=pedDat$X
	Z=pedDat$Z
		
	#if F1 TRUE, using identify working correlation
	#pw perturbation weight (Norm or Rade)
	Z<-scale(Z,scale=F) #Z need to be centered !
	#
	y<-as.vector(y)
	names(ID)=c("FID","IID","PAT","MAT")
	##sort according to FID
	p=ncol(Z); q=ncol(XC)-1; nx=q+1 
	#
	null<-gskat_fit(y,XC,ID,F1=F1)
	alpha<-null[[1]]; d<-null[[2]]
	mu<-as.vector(plogis(XC%*%alpha))      #mu<-as.vector(plogis(X%*%alpha))
	
	mydata<-data.frame(ID, X=XC, Z=Z, y=y, mu=mu)	
	result<-lapply(split(mydata,mydata$FID),function (adata) {
				yi<-adata[,(ncol(adata)-1)]; mui<-adata[,ncol(adata)]; Xi<-as.matrix(adata[,5:(5+nx-1)])	
				Zi<-as.matrix(adata[,(ncol(adata)-p-1):(ncol(adata)-2)])
				XZ<-cbind(Xi,Zi)
				
				if(length(mui)>1){Di<-diag((mui)*(1-mui))%*%XZ; Ai<-diag((mui)*(1-mui))} else {
					Di<-(mui)*(1-mui)*XZ; Ai<-(mui)*(1-mui)
				}	
				
				Rd<-mykin(adata$IID,adata$PAT,adata$MAT,F1=F1)
				if(length(Rd)>1){Rd<-Rd*d;diag(Rd)=1}
				Vi<-sqrt(Ai)%*%(Rd)%*%sqrt(Ai)
				iVi<-solve(Vi)
				
				Ui<-t(Zi)%*%Ai%*%iVi%*%t(t((yi-mui)))
				
				Covyi<-tcrossprod(yi-mui)
				#Delta%*%(X,Z)
				Bi<-t(Di)%*%iVi%*%Covyi%*%iVi%*%Di
				AAi<-t(Di)%*%iVi%*%Di
				return(list(Ui,Bi,AAi,Rd))
			})
	U=result[[1]][[1]];B=result[[1]][[2]];A=result[[1]][[3]]
	n=length(unique(ID$FID))
	for (i in 2:n) {
		U=U+result[[i]][[1]]
		B=B+result[[i]][[2]]
		A=A+result[[i]][[3]]
	}
	
	Azx<-A[(q+1+1):ncol(A),1:(q+1)]
	Axx<-A[1:(q+1),1:(q+1)]
	C<-cbind(-Azx%*%solve(Axx),diag(p))
	
	Ts<-t(U)%*%U
	
	#perturbation
	#R<-as.matrix(Matrix::bdiag(lapply(split(ID,ID$FID),function(adata){mykin(adata$IID,adata$PAT,adata$MAT,F1=F1)})))
	R<-as.matrix(bdiag(lapply(result,function(alist){return(alist[[4]])}))) #mom R
	a=diag((mu)*(1-mu))
	iV<-solve(sqrt(a)%*%(R)%*%sqrt(a))	
	#U_mat=t(Z)%*%a%*%iV%*%t(t((y-mu)))
	
	y_mu_t=t(t(y-mu)) 
	ZaV<-t(Z)%*%a%*%iV
	ni<-do.call(c,lapply(split(ID$FID,ID$FID),function(x){length(x)}))
	Ub.boot<-matrix(ncol=np,nrow=nrow(U))
	
	if (pw=="Norm") {
		Ub.boot<-apply(Ub.boot,2, function (x) {
					res.pert<-y_mu_t* rep(rnorm(n), ni) #Normal perturbation
					return(ZaV%*%res.pert)
				})		
	} else if (pw=="Rade") {
		Ub.boot<-apply(Ub.boot,2, function (x) {
					res.pert<-y_mu_t* rep(sample(c(1,-1),n,replace=T), ni)  #Rademacher distribution
					return(ZaV%*%res.pert)
				})		
	}
	
	Ub.boot<-t(Ub.boot) #
	Ts_boot<-apply(Ub.boot,1,function (x) {t(x)%*%x })
	
	mean(Ts_boot)
	df=12/kurtosis(Ts_boot)
	
	mu_Ts<- sum(diag(B%*%t(C)%*%C)) #trace(BC)  # Bsqrt=mysqrt(B); BC<-Bsqrt%*%t(C)%*%C%*%Bsqrt
	pval=1-pchisq((Ts-mu_Ts)*sqrt(2*df)/sqrt(var(Ts_boot))+df, df=df)
	var_Ts=2*(sum(diag(B%*%t(C)%*%C%*%B%*%t(C)%*%C)))
	return(list(pval=pval,Ts=Ts,mu_Ts=mu_Ts,var_Ts=var_Ts, PM=mean(Ts_boot),PV=var(Ts_boot)))
}


#########################################################################
#Auxiliary functions
########################################################################

mykin<-function(IID,PID,MID,F1=FALSE){
	#assume all ID characters
	#F1: if TRUE, forced to be identity matrix
	if (F1==TRUE){KIN2=diag(length(IID))} else if
	(length(IID)==1) {KIN2=1} else {	
		#require("kinship")
		ID=unique(c(IID,PID,MID)); ID=ID[ID!="0"]
		UID=ID[!ID%in%IID]
		IID<-c(IID,UID); PID=c(PID,rep(0,length(UID))); MID<-c(MID,rep(0,length(UID)))
		KIN2<-2*kinship(IID,PID,MID)[1:(length(IID)-length(UID)),1:(length(IID)-length(UID))]
		if(prod(diag(KIN2==1))!=1) {warning("Kinship is not correct")} 
	}
	return(KIN2)
}

blockMatrixDiagonal<-function(...){  
	matrixList<-list(...)
	if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]
	
	dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
	finalDimension<-sum(dimensions)
	finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
	index<-1
	for(k in 1:length(dimensions)){
		finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
		index<-index+dimensions[k]
	}
	finalMatrix
} 

mysqrt<-function (V) {
	#V=B
	V.eig<-eigen(V)
	V.eigvalues<-V.eig$values
	V.eigvalues[-1e-8<V.eigvalues & V.eigvalues<0]=0
	V.sqrt=V.eig$vectors %*% diag(sqrt(V.eigvalues)) %*% solve(V.eig$vectors)
	return(V.sqrt)
}


