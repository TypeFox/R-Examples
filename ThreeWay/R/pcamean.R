pcamean <-
function(X,n,m,p,laba,labb,labc){			

narg=nargs()
if (narg<5){
	laba=paste("a",1:n,sep="")
	labb=paste("b",1:m,sep="")
	labc=paste("c",1:p,sep="")
}

X=as.matrix(X)
check=1
while (check==1){
	cat("Over which mode do you want to compute means? Enter 1,2 or 3 (for A,B, or C)?",fill=TRUE)
	aggregmode=scan("",n=1)
	# 	 aggregmode = mode over which means are computed
	if (length(aggregmode)==0){
		aggregmode=0
	}
	if (aggregmode==1){
		check=0
		y1=m
		y2=p				# B-mode x C-mode
		check1=1
		cat(paste("A", m,"by", p," matrix for modes B x C will be computed"),fill=TRUE)
		while (check1==1)
			{
			cat("Which mode will serve as 'variables mode'? Enter 2 (=B) or 3 (=C):",fill=TRUE)
			loadingsmode=scan("",n=1)
			# loadingsmode = mode which is to be considered as the one 
			#  				for which you want to obtain simple structured loadings
			# Note that the other mode then is the "scores" mode 
			# with scores normalized to unit sums of squares
			if (length(loadingsmode)==0){
				loadingsmode=0
			}
			if (loadingsmode==2){ 
				check1=0
				load='P'
				cat("B will contain loadings, C will contain scores",fill=TRUE)
			}
			if (loadingsmode==3){
				check1=0
				load='Q'
				cat("C will contain loadings, B will contain scores",fill=TRUE)
			} else{
				if ((loadingsmode!=2) & (loadingsmode!=3)){
					cat(" ",fill=TRUE)
					cat("Error! Select a proper number!",fill=TRUE)
					cat(" ",fill=TRUE)
				}
			}
		}
	}
	if (aggregmode==2){
		check=0
		y1=p
		y2=n			# C-mode x A-mode
		X=permnew(X,n,m,p)		
		check1=1
		while (check1==1){
			cat(paste("A", p," by", n," matrix for modes C x A will be computed"),fill=TRUE)
			cat("Which mode will serve as 'variables mode'? Enter 1 (=A) or 3 (=C):",fill=TRUE)
			loadingsmode=scan("",n=1)
			if (length(loadingsmode)==0){
				loadingsmode=0
			}
			if (loadingsmode==3){
				check1=0
				load='P'
				cat("C will contain loadings, A will contain scores",fill=TRUE)
			}
			if (loadingsmode==1){
				check1=0
				load='Q'
				cat("A will contain loadings, C will contain scores",fill=TRUE)
			} else{
				if ((loadingsmode!=1) & (loadingsmode!=3)){
					cat(" ",fill=TRUE)
					cat("Error! Select a proper number!",fill=TRUE)
					cat(" ",fill=TRUE)
				}
			}
		}
	}
	if (aggregmode==3){
		check=0
		X=permnew(X,n,m,p)		
		X=permnew(X,m,p,n)		
		y1=n
		y2=m				# A-mode x B-mode
		check1=1
		while (check1==1){
			cat(paste("A", n," by", m," matrix for modes A x B will be computed"),fill=TRUE)
			cat("Which mode will serve as 'variables mode'? Enter 1 (=A) or 2 (=B):",fill=TRUE)
			loadingsmode=scan("",n=1)
			if (length(loadingsmode)==0){
				loadingsmode=0
			}
			if (loadingsmode==1){
				check1=0
				load='P'
				cat("A will contain loadings, B will contain scores",fill=TRUE)
			}
			if (loadingsmode==2){
				check1=0
				load='Q'
				cat("B will contain loadings, A will contain scores",fill=TRUE)
			} else{
				if ((loadingsmode!=1) & (loadingsmode!=2)){
					cat(" ",fill=TRUE)
					cat("Error! Select a proper number!",fill=TRUE)
					cat(" ",fill=TRUE)
				}
			}
		}
	}
	if ((aggregmode!=1) & (aggregmode!=2) & (aggregmode!=3)){
		cat(" ",fill=TRUE)
		cat("Error! Select a proper number!",fill=TRUE)
		cat(" ",fill=TRUE)
	}
}

mm=t(SUM(X)$mc)
Y=matrix(mm,y1,y2)  		# P-mode x Q-mode
if (nrow(Y)==n){
	rownames(Y)=laba
}
if (ncol(Y)==n){
	colnames(Y)=laba
}
if (nrow(Y)==m){
	rownames(Y)=labb
}
if (ncol(Y)==m){
	colnames(Y)=labb
}
if (nrow(Y)==p){
	rownames(Y)=labc
}
if (ncol(Y)==p){
	colnames(Y)=labc
}
if (sum(Y^2)<1e-10){
	cat("Means to analyze are exteremely small, probably data have been centered across the",fill=TRUE)
	cat("same mode as over which means are computed, which makes no sense.",fill=TRUE)
	cat("Therefore, PCA of means is not carried out!",fill=TRUE)
	A1=NULL
	A2=NULL
	B1=NULL
	B2=NULL
	C1=NULL
	C2=NULL
    ev=NULL
} else{
	cat("If you want to standardize 'variables', type '1':",fill=TRUE)
	cc=scan("",n=1)
	if (length(cc)==0){
		cc=0
	}
	if (cc==1){
		if (load=='Q'){
			Y=nrm2(Cc(Y))
		}
		if (load=='P'){
			Y=t(nrm2(Cc(t(Y))))
		}
		cat("'Variables' have been standardized",fill=TRUE)
	} else{
		cat("'Variables' have not been standardized",fill=TRUE)
	}
}

if (sum(Y^2)>1e-10){
SVD=svd(Y)		 # P for rows  Q for columns
ev=SVD$d^2
rr=length(ev)
cat("Eigenvalues and cumulative percentages (Scree plot in Figure)",fill=TRUE)
mat=matrix(,length(ev),2)		
mat[,1]=ev
mat[,2]=cumsum(ev)/sum(ev)*100
labComp=paste("Comp.",1:rr,sep="")
names(ev)=labComp
rownames(mat)=labComp
colnames(mat)=c("Eigenvalue","Fit(%)")
print(round(mat,digits=2))
plot(t(c(1:rr)),ev,type="l",xlab="",ylab="")

rMax=min(nrow(Y),ncol(Y))
while(check==0){
	cat("How many components do you want to maintain for rotation?",fill=TRUE)
	r=scan("",n=1)
	# check on the number of components
	if (length(r)==0){
		r=0
	}
	if ((r>1) & (r<=rMax) & ((floor(r)-r)==0)){
		check=1
		labComp=paste("Comp.",1:r,sep="")
		if (load=='P'){
			# varimax of P-mode loadings
			NRMVAR=normvari(SVD$u[,1:r]%*%diag(SVD$d,nrow=rMax)[1:r,1:r])			
			Q1=SVD$v[,1:r]%*%NRMVAR$T
			sg=diag(sign(colSums(NRMVAR$B)),nrow=ncol(NRMVAR$B))			# reflect if necessary
			P1=NRMVAR$B%*%sg
			Q1=Q1%*%sg
			# Harris-Kaiser of P-mode loadings
			NRMVAR2=normvari(SVD$u[,1:r])					
			Q2=SVD$v[,1:r]%*%diag(SVD$d,nrow=rMax)[1:r,1:r]%*%NRMVAR2$T
			Ds=diag(colSums(Q2^2),nrow=ncol(Q2))
			P2=NRMVAR2$B%*%Ds^.5
			Q2=Q2%*%solve(Ds^.5)		
			sg=diag(sign(colSums(P2)),nrow=ncol(P2))			 # reflect if necessary
			P2=P2%*%sg
			Q2=Q2%*%sg
			cat("A1, B1 and/or C1 contain solution based on varimax rotation of loadings",fill=TRUE)
			cat("A2, B2 and/or C2 contain solution based on oblique 'HKIC' orthomax rotation of loadings",fill=TRUE)
		} 
		if (load=='Q'){
			# varimax of Q-mode loadings
			NRMVAR=normvari(SVD$v[,1:r]%*%diag(SVD$d,nrow=rMax)[1:r,1:r])
			P1=SVD$u[,1:r]%*%NRMVAR$T
			sg=diag(sign(colSums(NRMVAR$B)),nrow=ncol(NRMVAR$B))		# reflect if necessary
			P1=P1%*%sg
			Q1=NRMVAR$B%*%sg
			# Harris-Kaiser of Q-mode loadings
			NRMVAR2=normvari(SVD$v[,1:r])
			P2=SVD$u[,1:r]%*%diag(SVD$d,nrow=rMax)[1:r,1:r]%*%NRMVAR2$T
			Ds=diag(colSums(P2^2),nrow=ncol(P2))
			Q2=NRMVAR2$B%*%Ds^.5
			P2=P2%*%solve(Ds^.5)
			sg=diag(sign(colSums(Q2)),nrow=ncol(Q2))			 # reflect if necessary
			P2=P2%*%sg
			Q2=Q2%*%sg
			cat("A1, B1 and/or C1 contain solution based on varimax rotation of loadings",fill=TRUE)
			cat("A2, B2 and/or C2 contain solution based on oblique 'HKIC' orthomax rotation of loadings",fill=TRUE)
		}
	} else{
		if (r==1){
			labComp=c("Comp.1")
			check=1
			if (load=='P'){
				P1=SVD$u[,1]*SVD$d[1]
				sg=sign(sum(P1))			# reflect if necessary
				P1=as.matrix(P1)*sg
				Q1=as.matrix(SVD$v[,1])*sg
				P2=P1
				Q2=Q1
				cat("A1, B1 and/or C1 contain principal component solution",fill=TRUE)
			} 
			if (load=='Q'){
				Q1=SVD$v[,1]*SVD$d[1]
				sg=sign(sum(Q1))		# reflect if necessary
				Q1=as.matrix(Q1)*sg
				P1=as.matrix(SVD$u[,1])*sg
				P2=P1
				Q2=Q1
				cat("A1=A2, B1=B2 and/or C1=C2 contain principal component solution",fill=TRUE)
			}
		} else{
			cat(" ",fill=TRUE)
			cat("Error! Select a proper number!",fill=TRUE)
			cat(" ",fill=TRUE)
		}
	}
}     

if (aggregmode==1 ){
	A1=NULL
	A2=NULL
	B1=P1
	B2=P2
	C1=Q1
	C2=Q2
}
if (aggregmode==2){
	A1=Q1
	A2=Q2
	B1=NULL
	B2=NULL
	C1=P1
	C2=P2
}
if (aggregmode==3){
	A1=P1
	A2=P2
	B1=Q1
	B2=Q2
	C1=NULL
	C2=NULL
}
if (identical(A1,NULL)==FALSE){
	rownames(A1)=laba
	rownames(A2)=laba
	colnames(A1)=labComp
	colnames(A2)=labComp
}
if (identical(B1,NULL)==FALSE){
	rownames(B1)=labb
	rownames(B2)=labb
	colnames(B1)=labComp
	colnames(B2)=labComp
}
if (identical(C1,NULL)==FALSE){
	rownames(C1)=labc
	rownames(C2)=labc
	colnames(C1)=labComp
	colnames(C2)=labComp
}
}

out=list()
out$Y=Y
out$ev=ev
out$A1=A1
out$A2=A2
out$B1=B1
out$B2=B2
out$C1=C1
out$C2=C2	
return(out)
}
