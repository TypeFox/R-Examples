# TODO: GSKAT main function block (rare variants)
# Beta version!
# Author: Xuefeng Wang
###############################################################################


#################################BURDEN TEST##################################
genoC<-function(Z,w_a=1,w_b=5,BM="WST") {
	maf<-apply(Z,2,sum)/(nrow(Z)*2)
	#if (model=="R"){Z[Z==1]=0}  #recessive  	
	
	if (BM=="WST"){
		#wi<- (dbeta(maf,w_a, w_b))^2  #beta density
		wi<- dbeta(maf,w_a,w_b)
		Z1<- Z%*%t(t(wi))   #combined gentoypes
	} else if (BM=="CAST"){
		Zr<-Z[,maf<0.03]     #rare, threshold 0.03
		if (ncol(Zr)<1){stop("no rare")}
		Z1<-rowSums(Zr)
		Z1<-as.numeric(Z1>0)
	} else if (BM=="count"){
		Zr<-Z[,maf<0.03]     #rare, threshold 0.03
		if (ncol(Zr)<1){stop("no rare")}
		Z1<-apply(Zr,1,function (arow) {
					sum(arow>0)
				})
	}
	return(Z1)
}
	

##null fit function: estimate aplha.hat and delta.hat
myfit.RFAM_IC<-function(y,XC,FID){
	#y:phenotype
	#XC:covaraite matrix, including PC and intercept
	#FID: family ID
	
	#alpha<-t(t(c(-1,rep(1,(ncol(XC)-1))))) #initial values
	alpha<-summary(glm(y~XC[,-1],family=binomial("logit")))$coeff[,1]
	#OR use glm to get initial values	
	#if (class(FID)=="factor"){FID=as.character(FID)}
	#n<-length(unique(FID)) # family number
	#R1<-matrix(0.5,ncol=4,nrow=4); diag(R1)=1; R1[1,2]=R1[2,1]=0 # 2*Kinship matrix
	#R1<-diag(4)
	
	diff1<-10
	i<-0
	while (diff1>10e-6){
		
		i<-i+1
		
		#Fisher scoring
		mu<- plogis( XC %*% alpha)
		mu=as.vector(mu)
		#D=diag((mu)*(1-mu))%*%X
		V=diag((mu)*(1-mu)) #V=sqrt(A)%*%(R)%*%sqrt(A) 
		U=t(XC)%*%t(t((y-mu)))    #U=t(D)%*%solve(V)%*%t(t((y-mu)))
		I=t(XC)%*%V%*%XC	 	#I<-t(D)%*%solve(V)%*%D
		
		alpha.new<-alpha+solve(I)%*%U
		diff1 <- max(abs(alpha.new-alpha))
		alpha<-alpha.new
		
#		#Method of moments
#		mu.new<-plogis(X%*%alpha)
#		r.hat<-(y-mu.new)/sqrt((mu.new*(1-mu.new)))    #Pearson residual
#		
#		mydata<-data.frame(FID=rep(1:n,each=4),r.hat=r.hat)
#		
#		result<-Reduce("+",lapply(split(mydata,mydata$FID),function (adata) {
#			ri<-adata[,2]
#			RI<-tcrossprod(ri)
#			mynu<-sum(RI[upper.tri(RI)])
#			myde<-sum(R1[upper.tri(R1)])
#			return(c(mynu,myde))
#		}))
#		d.new<-result[1]/result[2]
#		diff2 <- abs(d.new - d)		
#		d<-d.new
	}	
	return(list(alpha,i,mu,V))
}

##score test
score.RFAM_IC_burden<-function(y,XC,Z,FID) {
	#y: binary phenotype coded as 0, 1
	#XC:covaraite matrix, including PC and intercept 
	#Z: combined genotype matrix (one col.)
	
	#require("CompQuadForm")	
	#
	XC=cbind(1,data.matrix(XC)) #now with intercept
	Z<-data.matrix(Z)	#Z is the combined genotype, must be matrix format
	y<-as.vector(y);FID<-as.factor(FID)
	##sort according to FID	
	FID.old<-FID
	FID<-FID[order(FID.old)]; y<-y[order(FID.old)];XC<-XC[order(FID.old),];Z<-data.matrix(Z[order(FID.old),])
	p=ncol(Z); q=ncol(XC)-1
	#
	null<-myfit.RFAM_IC(y,XC,FID)
	alpha<-null[[1]]
	#
	#mu<-as.vector(null[[3]]) 
	mu<-as.vector(plogis(XC%*%alpha))
	U=t(Z)%*%t(t((y-mu)))
	#V=null[[4]] 
	V=diag((mu)*(1-mu))
	Covy<-tcrossprod(y-mu)
	#Covy[((row(Covy)-1)%/%4+1)!=((col(Covy)-1)%/%4+1)]=0  #for same family strut of 4 members
	Covy=Covy*blockMatrixDiagonal(lapply(split(FID,FID),function(alist){matrix(1,length(alist),length(alist))}))
	
	XZ<-cbind(XC,Z)
	B<-t(XZ)%*%Covy%*%XZ  #Bi<-t(Di)%*%iVi%*%Covyi%*%iVi%*%Di
	A<-t(XZ)%*%V%*%XZ  #AAi<-t(Di)%*%iVi%*%Di 
	
	Azx<-A[(q+1+1):ncol(A),1:(q+1)]
	Axx<-A[1:(q+1),1:(q+1)]
	C<-cbind(-Azx%*%solve(Axx),diag(p))
	#foo<<-p
	
	TS<-t(U)%*%U
	Bsqrt=mysqrt(B)
	BC<-Bsqrt%*%t(C)%*%C%*%Bsqrt
	
	Lamda<-eigen(BC,only.values=T)$values
	davies(TS,Lamda,rep(1, length(Lamda)))$Qq
}


gee_wrap<-function(y,XC,Z,FID,corstr="independence"){
	#XC:covaraites and PC, no intercept
	#Z: one col matrix
	#require(gee)
	y<-data.matrix(y);XC<-data.matrix(XC);Z<-data.matrix(Z)
	FID<-as.factor(FID)
	FID.old<-FID
	FID<-FID[order(FID.old)]; y<-y[order(FID.old),];XC<-XC[order(FID.old),];Z<-Z[order(FID.old),] #!!
	
	gfit= try(gee(y~XC+Z,id=FID,family=binomial,corstr=corstr,silent=T))
	if (class(gfit)[1]=="try-error"){pv2<-9999} else {
		coeff=summary(gfit)$coefficients
		z_score=coeff[nrow(coeff),ncol(coeff)] #z-score of Z
		pv2=2*(pnorm(abs(z_score),lower.tail=F))
	}
	pv2
}


geeglm_wrap<-function(y,XC,Z,FID,corstr="independence"){
	#y:phenotyep matrix,
	#XC:covaraites and PC, no intercept
	#Z: one col matrix
	#FID: need to be ordered in the following
	#require(geepack)
	y<-data.matrix(y);XC<-data.matrix(XC);Z<-data.matrix(Z)
	FID<-as.factor(FID)
	FID.old<-FID
	FID<-FID[order(FID.old)]; y<-y[order(FID.old),];XC<-XC[order(FID.old),];Z<-Z[order(FID.old),] #!!
	fit=geeglm(y~XC+Z,id=FID,family=binomial,corstr=corstr)
	coeff=summary(fit)$coefficients # estimates and p-vals
	pout<-coeff[nrow(coeff),ncol(coeff)] 
	pout
}
#geeglm_wrap(y,XC=X,Z,FID)


#####################################GEE_SKAT###############################
##null fit function: estimate aplha.hat and delta.hat
fit_FSKAT_IC<-function(y,XC,FID){
	#y:phenotype
	#XC:covaraite matrix, including PC and intercept !
	#FID: family ID, not used if indentity working correlation
	
	#alpha<-t(t(c(-1,rep(1,(ncol(XC)-1))))) #initial values
	alpha<-summary(glm(y~XC[,-1],family=binomial("logit")))$coeff[,1]
	#OR use glm to get initial values	
	#if (class(FID)=="factor"){FID=as.character(FID)}
	#n<-length(unique(FID)) # family number
	#R1<-matrix(0.5,ncol=4,nrow=4); diag(R1)=1; R1[1,2]=R1[2,1]=0 # 2*Kinship matrix
	#R1<-diag(4)
	
	diff1<-10
	i<-0
	while (diff1>10e-6){
		
		i<-i+1
		
		#Fisher scoring
		mu<- plogis( XC %*% alpha)
		mu=as.vector(mu)
		#D=diag((mu)*(1-mu))%*%X
		V=diag((mu)*(1-mu)) #V=sqrt(A)%*%(R)%*%sqrt(A) 
		U=t(XC)%*%t(t((y-mu)))    #U=t(D)%*%solve(V)%*%t(t((y-mu)))
		I=t(XC)%*%V%*%XC	 	#I<-t(D)%*%solve(V)%*%D
		
		alpha.new<-alpha+solve(I)%*%U
		diff1 <- max(abs(alpha.new-alpha))
		alpha<-alpha.new
		
#		#Method of moments
#		mu.new<-plogis(X%*%alpha)
#		r.hat<-(y-mu.new)/sqrt((mu.new*(1-mu.new)))    #Pearson residual
#		
#		mydata<-data.frame(FID=rep(1:n,each=4),r.hat=r.hat)
#		
#		result<-Reduce("+",lapply(split(mydata,mydata$FID),function (adata) {
#			ri<-adata[,2]
#			RI<-tcrossprod(ri)
#			mynu<-sum(RI[upper.tri(RI)])
#			myde<-sum(R1[upper.tri(R1)])
#			return(c(mynu,myde))
#		}))
#		d.new<-result[1]/result[2]
#		diff2 <- abs(d.new - d)		
#		d<-d.new
	}	
	return(list(alpha,i))
}

##score test
score_FSKAT_IC_pertu<-function(y,XC,Z,ID,w_a=1,w_b=5,pw="Rade",Uc=FALSE,sW=TRUE,np=10000) {
	#w_a=1;w_b=25; pw="Rade"; Uc=FALSE
	maf<-apply(Z,2,sum)/(nrow(Z)*2) #cacluate before center
	Z<-scale(Z,scale=F) #Z need to be centered !!!
	#require("CompQuadForm")	
	#require(e1071)
	#
	XC=cbind(1,data.matrix(XC)) #now with intercept
	Z<-data.matrix(Z)	#Z is the combined genotype, must be matrix format
	y<-as.vector(y);FID<-as.factor(ID$FID)
	##sort according to FID	
	FID.old<-FID
	FID<-FID[order(FID.old)]; y<-y[order(FID.old)];XC<-XC[order(FID.old),];Z<-data.matrix(Z[order(FID.old),]) #!!
	ID<-ID[order(FID.old),]
	p=ncol(Z); q=ncol(XC)-1
	#
	
	#weight W
	#w_a=1; w_b=25
	if (length(maf)==1){W=1} else {
		W<-diag((dbeta(maf,w_a, w_b))^2) } #beta density
	if (sW==TRUE) {W<-W/max(W)}
	#W<-diag(length(maf))
	
	#
	null<-fit_FSKAT_IC(y,XC,FID)
	alpha<-null[[1]]
	#
	mu<-as.vector(plogis(XC%*%alpha))
	U=t(Z)%*%t(t((y-mu)))
	#V=null[[4]] 
	V=diag((mu)*(1-mu))
	Covy<-tcrossprod(y-mu)
	#Covy[((row(Covy)-1)%/%4+1)!=((col(Covy)-1)%/%4+1)]=0  #for same family strut of 4 members
	Covy=Covy*blockMatrixDiagonal(lapply(split(FID,FID),function(alist){matrix(1,length(alist),length(alist))}))
	
	XZ<-cbind(XC,Z)
	B<-t(XZ)%*%Covy%*%XZ  #Bi<-t(Di)%*%iVi%*%Covyi%*%iVi%*%Di
	A<-t(XZ)%*%V%*%XZ  #AAi<-t(Di)%*%iVi%*%Di 
	
	Azx<-A[(q+1+1):ncol(A),1:(q+1)]
	Axx<-A[1:(q+1),1:(q+1)]
	C<-cbind(-Azx%*%solve(Axx),diag(p))
	
	TS<-t(U)%*%W%*%U
	
	Bsqrt=mysqrt(B)
	BC<-Bsqrt%*%t(C)%*%W%*%C%*%Bsqrt	
	Lamda<-eigen(BC,only.values=T)$values
	results<-davies(TS,Lamda,rep(1, length(Lamda)))
	#results
	
	#perturbation
	#mom R or Identity
	a=diag((mu)*(1-mu))
	R<-diag(dim(a)[1])
	iV<-solve(sqrt(a)%*%(R)%*%sqrt(a))	
	#U_mat=t(Z)%*%a%*%iV%*%t(t((y-mu)))
	
	y_mu_t=t(t(y-mu)) 
	ZaV<-t(Z)%*%a%*%iV
	ni<-do.call(c,lapply(split(ID$FID,ID$FID),function(x){length(x)}))
	Ub.boot<-matrix(ncol=np,nrow=nrow(U))
	n=length(unique(ID$FID))
	
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
	
	if (Uc==TRUE){Ub.boot<-apply(Ub.boot,1,scale, scale = FALSE)} else if (length(maf)==1) {
		Ub.boot<-t(t(Ub.boot))} else {
		Ub.boot<-t(Ub.boot)} #
	Ts_boot<-apply(Ub.boot,1,function (x) {t(x)%*%W%*%x })
	
	mean(Ts_boot)
	var(Ts_boot); var_theory=2*(sum(diag(B%*%t(C)%*%W%*%C%*%B%*%t(C)%*%W%*%C)))
	df=12/kurtosis(Ts_boot)
	if (df<0) {df=100}
	
	mu_Ts<- sum(diag(B%*%t(C)%*%W%*%C)) #trace(BC)  # Bsqrt=mysqrt(B); BC<-Bsqrt%*%t(C)%*%C%*%Bsqrt
	#mu_Ts
	if (df<0) {pval=NA} else {
		pval=1-pchisq((TS-mu_Ts)*sqrt(2*df)/sqrt(var(Ts_boot))+df, df=df)}
	#pval
	
	return(list(pval_davies=results$Qq,ifault=results$ifault,PM=mean(Ts_boot),mu_Ts=mu_Ts, PV=var(Ts_boot),var_Ts=var_theory ,pval_pert=pval))
}









