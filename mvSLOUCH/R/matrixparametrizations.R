## ------------- parametrization of symmetric positive definite matrix ------------------------------
.sym.par <- function (x,nchar) {
## function taken from ouch package
    mSym<-NA
    nchar <- floor(sqrt(2*length(x)))
    if (nchar*(nchar+1)!=2*length(x)) {
        print("a symmetric matrix is parameterized by a triangular number of parameters")
    }
    y <- matrix(0,nchar,nchar)
    y[lower.tri(y,diag=TRUE)] <- x
    mSym<-y%*%t(y)
    mSym
}

.sym.unpar <- function (x) {
## function taken from ouch package
  y <- t(chol(x))
  y[lower.tri(y,diag=TRUE)]
}

.par.transform.symmetric.matrix<-function(vParams,n){
    A<-matrix(NA,n,n)
    A[upper.tri(A,diag=TRUE)]<-vParams
    A[lower.tri(A)]<-t(A)[lower.tri(A)]
    A
}

.par.inv.transform.symmetric.matrix<-function(A){
    A[upper.tri(A,diag=TRUE)]
}
## -------------------------------------------------------------------------------------------------





## ------------- Givens rotation functions ---------------------------------------------------------
## these functions are taken nearly as is from :
## for the real case :
## G. Golub, C. Van Loan; Matrix Computations; The Johns Hopkins University Press; 1989; Baltimore
## for the complex case : ## this still has to be done and worked on
## D. Bindel, J. Demmel, W. Kahan; On Computing Givens Rotations Reliably and Efficiently;ACM Transactions on Mathematical Software 28(2):206-238

.givens.real<-function(a,b,getrho=TRUE,cs=NULL){## algorithm 5.1.5 p202 and eq 5.1.9 p204
    if (is.null(cs)){
	if (b==0){cG<-1;sG<-0}
	else{
	    if (abs(b)>abs(a)){tau<-(-1)*a/b;sG<-1/sqrt(1+tau^2);cG<-sG*tau}
	    else{tau<-(-1)*b/a;cG<-1/sqrt(1+tau^2);sG<-cG*tau}
	}
    }else{cG<-cs[1];sG<-cs[2];}
    if (getrho){
	if (cG==0){rho<-1}
	else{
	    if(abs(sG)<abs(cG)){rho<-sign(cG)*sG/2}
	    else{rho<-2*sign(sG)/cG}	    
	}	
    }else{rho<-c(cG,sG)}
    rho
}

.inv.rho.givens.real<-function(rho){## eq 5.1.10 p205 -> to get a and b back do .row.rot or .col.rot
    if (rho==1){cG<-0;sG<-1}
    else{
	if (abs(rho)<1){sG<-2*rho;cG<-sqrt(1-sG^2)}
	else{cG<-2/rho;sG<-sqrt(1-cG^2)}
    }
    c(cG,sG)
}

.row.rot<-function(M,cG,sG){## algorithm 5.1.6 p203
    if (nrow(M)!=2){print("Cannot perform .row.rot when M does not have 2 rows");rotM<-NA}
    else{rotM<-apply(M,2,function(x){c(cG*x[1]-sG*x[2],sG*x[1]+cG*x[2])})}
    rotM
}

.col.rot<-function(M,cG,sG){## algorithm 5.1.7 p203
    if (ncol(M)!=2){print("Cannot perform .col.rot when M does not have 2 columns");rotM<-NA}
    else{rotM<-apply(M,1,function(x){c(cG*x[1]-sG*x[2],sG*x[1]+cG*x[2])})}
    rotM
}

.givens.mult.M<-function(M,cG,sG,i1,i2){
    M[c(i1,i2),]<-.row.rot(M[c(i1,i2),],cG,sG)
    M
}

.M.mult.givens<-function(M,cG,sG,j1,j2){
    M[,c(j1,j2)]<-.col.rot(M[,c(j1,j2)],cG,sG)
    M
}
## -------------------------------------------------------------------------------------------------

## ------------ parametrization of traingular matrix -----------------------------------------------
.par.transform.uppertri.matrix<-function(vParams,k){
    A<-matrix(0,k,k)
    A[upper.tri(A,diag=TRUE)]<-vParams
    A
}

.par.inv.transform.uppertri.matrix<-function(A){
    A[upper.tri(A,diag=TRUE)]
}

.par.transform.lowertri.matrix<-function(vParams,k){
    A<-matrix(0,k,k)
    A[lower.tri(A,diag=TRUE)]<-vParams
    A
}

.par.inv.transform.lowertri.matrix<-function(A){
    A[lower.tri(A,diag=TRUE)]
}
## -------------------------------------------------------------------------------------------------



## ------------ parametrization of invertible matrix -----------------------------------------------
## ------------ needed for A's 2x2  parametrization ------------------------------------------------
.par.transform.invert.matrix<-function(vParams,n,method="qr",minRdiag=0.01){
    lA<-NA
    if (method=="qr"){lA<-.par.transform.invert.qr.matrix(vParams,n,NULL,NULL,minRdiag)}
    if (method=="svd"){lA<-.par.transform.invert.svd.matrix(vParams,n)}
    lA
}

.par.inv.transform.invert.matrix<-function(A,tol=1e-15,method="qr",mCsigns=NULL,minRdiag=0.01){
    lAparams<-NA
    if (method=="qr"){lAparams<-.par.inv.transform.invert.qr.matrix(A,mCsigns,tol,minRdiag)}
    if (method=="svd"){lAparams<-.par.inv.transform.invert.svd.matrix(A,tol)}
    lAparams
}

.par.transform.invert.svd.matrix<-function(vParams,n){
## parametrization is based on SVD but not exactly : D is allowed to be negative
## the reason is for the parametrization of U,V on the diagonal after Givens zeroing
## we can have negative values, ie rho does not hold information on the sign
    A<-NA
    if (length(vParams)!=n^2){print(paste("Cannot calculate A, vParams should be of length : ",n^2))}
    else{
	k<-(length(vParams)-n)/2
	D<-diag(vParams[c(1:n)],n,n)
	U<-.par.transform.orth.matrix.givens(vParams[(n+1):(n+k)],n)
	V<-.par.transform.orth.matrix.givens(vParams[(n+k+1):(n^2)],n)
	A<-U%*%D%*%t(V)	
    }
    A
}

.par.inv.transform.invert.svd.matrix<-function(A,tol=1e-15){
    n<-nrow(A)
    vParams<-rep(NA,n^2)
    svdA<-svd(A)    
    vParamsU<-.par.inv.transform.orth.matrix.givens(svdA$u,tol,TRUE)
    vParamsV<-.par.inv.transform.orth.matrix.givens(svdA$v,tol,TRUE)
    vParams[1:n]<-svdA$d*vParamsU[1:n]*vParamsV[1:n] ## this allows the ``singular values'' to be negative but we don't care
    vParamsU<-vParamsU[-(1:n)]
    vParamsV<-vParamsV[-(1:n)]
    vParams[(n+1):(n^2)]<-c(vParamsU,vParamsV)
    vParams    
}

.par.transform.twobytwo.matrix<-function(vParams){
    matrix(vParams,2,2)
}

.par.inv.transform.twobytwo.matrix<-function(A){
    c(A)
}

.par.transform.invert.qr.matrix<-function(vParams,n,mCsigns=NULL,vDiagQ=NULL,minRdiag=0.01){
## the parametrization is done via QR decomposition assuming Rs diagonal is meant to be positive
    A<-NA
    mQcSigns<-NULL
    if(!is.null(mCsigns)){mQcSigns<-mCsigns}
    if (length(vParams)!=n^2){print(paste("Number of parameters describing A must be ",n,"^2=",n^2))}
    else{
	if (n==1){A<-vParams[1]}
	else{
	    k<-n*(n-1)/2
	    ## first n(n-1)/2 parameters are the Givens rotions of Q
	    lQparam<-.par.transform.orth.matrix.givens(vParams[1:k],n,mQcSigns,vDiagQ)
	    Q<-lQparam$Q
	    mQcSigns<-lQparam$cSigns
	    ## next n parameters are R's diagonal, they must be positive
	    R<-matrix(0,n,n)
	    diag(R)<-exp(vParams[(k+1):(k+n)])+minRdiag
	    ## last n(n-1)/2 parameters are R's upper triangle
	    R[upper.tri(R)]<-vParams[(k+n+1):(n^2)]	
	    A<-Q%*%R*(-1) ## turns out in R's qr need to multiply by (-1) to get back original
	}
    }
    list(A=A,QcSigns=mQcSigns)
}

.par.inv.transform.invert.qr.matrix<-function(A,mCSigns=NULL,tol=1e-15,minRdiag=0.01){
    n<-nrow(A)    
    vParams<-rep(NA,n^2)
    if (ncol(A)!=n){print("A must be a square matrix")}
    else{
	if (n==1){vParams[1]<-A[1,1];mCSigns=1}
	else{
    	    k<-n*(n-1)/2
	    qrA<-qr(A)   		
	    QqrA<-qr.Q(qrA)
	    RqrA<-qr.R(qrA)
    	    minRdiag<-min(minRdiag,min(diag(RqrA)/2))
	    vNegR<-which(diag(RqrA<0))
    	    RqrA[vNegR,]<-RqrA[vNegR,]*(-1)
	    QqrA[,vNegR]<-QqrA[,vNegR]*(-1)
	    lQparams<-.par.inv.transform.orth.matrix.givens(QqrA,mCSigns,tol)
	    vParams[1:k]<-lQparams$vParams
	    mCSigns<-lQparams$cSigns
	    vDiagQ<-lQparams$vDiagQ
	    vParams[(k+1):(k+n)]<-log(diag(RqrA)-minRdiag)
	    vParams[(k+n+1):(n^2)]<-RqrA[upper.tri(RqrA)]
	}	    
    }
    list(vParams=vParams,cSigns=mCSigns,vDiagQ=vDiagQ)
}
## -------------------------------------------------------------------------------------------------



## ------------ parametrization of decomposable matrix ---------------------------------------------
## ----------------- symmetric positive definite matrix --------------------------------------------
## this is meant to be an alternative parametrization to the one in sym.par (from ouch)

.par.transform.decomp.sym.matrix<-function(vParams,n){
## the eigenvalue decomposition of A=PJP^(T) P parametrized by givens rotations
## it is assumed that the eigenvectors (P's columns )are of unit length 
    A<-NA
    if (length(vParams)!=n^2){print(paste("Cannot calculate A, vParams should be of length : ",n^2))}
    else{
	if (n==1){A<-matrix(vParams[1],1,1)}
	else{
	    k<-n*(n-1)/2
	    if (is.complex(vParams[1:n])){
	    	print("Thera are complex eigenvalues, calling complex calculations!")
		A<-.par.transform.decomp.matrix(vParams,n)
	    }else{
		J<-diag(vParams[1:n],n,n) ## the eigenvalues	    
		P<-.par.transform.orth.matrix.givens(vParams[(n+1):(n+k)],n)
		A<-P%*%J%*%t(P)
	    }
	}
    }
    A
}

.par.inv.transform.decomp.sym.matrix<-function(A,tol=1e-15){
    n<-nrow(A)    
    vParams<-rep(NA,n^2)
    if (ncol(A)!=n){print("A must be a square matrix")}
    else{
	if (n==1){vParams<-A[1,1]}
	else{
	    eigA<-eigen(A)
	    if(is.complex(eigA$values)){
		print("A has complex eigen values, calling complex parametrization!")
		vParams<-.par.inv.transform.decomp.matrix(A)
	    }else{
	        k<-n*(n-1)/2
		vParams[1:n]<-eigA$values
		vParams[(n+1):(n+k)]<-.par.inv.transform.orth.matrix.givens(eigA$vectors,tol)		
	    }	    
	}
    }
    vParams
}


## ----------------- real eigenvalues --------------------------------------------------------------

## ---------------- positive definite --------------------------------------------------------------
.par.transform.decomp.pos.real.matrix<-function(vParams,n,mCSigns=NULL,minEigenValue=0.01){
## first n parameters are the logs of the eigenvalues -> we want all to be positive
## next numbers are somehow the eigenvectors
    vParams[1:n]<-exp(rev(sort(vParams[1:n])))+minEigenValue
    .par.transform.decomp.real.matrix(vParams,n,mCSigns)    
}

.par.inv.transform.decomp.pos.real.matrix<-function(A,eigenSigns=NULL,mCSigns=NULL,tol=1e-15,minEigenValue=0.01){
    lAparams<-.par.inv.transform.decomp.real.matrix(A,eigenSigns,mCSigns,tol)
    n<-nrow(A)
    lAparams$vParams[1:n]<-log(lAparams$vParams[1:n]-minEigenValue)
    lAparams
}

.par.transform.decomp.neg.real.matrix<-function(vParams,n,mCSigns=NULL,maxEigenValue=0.01){
## first n parameters are the logs of the eigenvalues -> we want all to be positive
## next numbers are somehow the eigenvectors
    vParams[1:n]<-(-1)*exp(rev(sort(vParams[1:n])))-maxEigenValue
    .par.transform.decomp.real.matrix(vParams,n,mCSigns)    
}

.par.inv.transform.decomp.neg.real.matrix<-function(A,eigenSigns=NULL,mCSigns=NULL,tol=1e-15,maxEigenValue=0.01){
    lAparams<-.par.inv.transform.decomp.real.matrix(A,eigenSigns,mCSigns,tol)
    n<-nrow(A)
    lAparams$vParams[1:n]<-log((-1)*lAparams$vParams[1:n]-maxEigenValue)
    lAparams
}

.par.transform.decomp.real.matrix<-function(vParams,n,mCSigns=NULL){
## the eigenvalue decomposition of A=PJP^(-1) P parametrized by QR 
## it is assumed that the eigenvectors (P's columns )are of unit length -> this should make
## R's columns also of unit length -> we assume diagonal of R is positive!
## here we do not care for the diagonal of Q after the Givens rotations as this will cancel out
## in PJP^-1
    A<-NA
    eigenSigns<-NA
    if (length(vParams)!=n^2){print(paste("Cannot calculate A, vParams should be of length : ",n^2))}
    else{
	if (n==1){A<-matrix(vParams[1],1,1);eigenSigns<-c(1)}
	else{	    
	    if (is.complex(vParams[1:n])){
	    	print("Thera are complex eigenvalues, calling complex calculations!")
		A<-.par.transform.decomp.matrix(vParams,n)
	    }else{
		k<-n*(n-1)/2
		J<-diag(vParams[1:n],n,n) ## the eigenvalues	    
		if (length(which(diag(J)==0))>0){print("WARNING : we have zero eigenvalues, matrix is singular")}
		lQparam<-.par.transform.orth.matrix.givens(vParams[(n+1):(n+k)],n,mCSigns)
		Q<-lQparam$Q
		mQcSigns<-lQparam$cSigns
		mR<-matrix(0,n,n)
		mR[upper.tri(mR)]<-vParams[(n+k+1):(n^2)]
		mR<-sapply(1:ncol(mR),function(i,mR){
		    k<-nrow(mR)
		    orgcolR<-mR[,i]
		    colR<-c(rep(NA,i-1),rep(0,n-i+1))
		    if (i>1){
		    colR[1]<-atan(orgcolR[1])*2/pi
		    if (i>2){
			for (j in 2:(i-1)){
			    a<-sqrt(1-sum(colR[1:(j-1)]^2))
			    colR[j]<-atan(orgcolR[j])*2*a/pi
			}
		    }
		    }
		    colR
		},mR=mR,simplify=TRUE)
		diag(mR)<-c(1,sapply(2:ncol(mR),function(j){sqrt(1-sum(mR[,j]^2))})) ## should work ok as we previously have zeroed R
		P<-Q%*%mR
		eigenSigns<-sign(P[1,])
		invP<-solve(mR)%*%t(Q) ## inversing using the QR would be quicker and faster but we need to fix the eigenvectors
		A<-P%*%J%*%invP	
	    }
	}
    }
    list(A=A,eigenSigns=eigenSigns,QcSigns=mQcSigns)
}

.par.inv.transform.decomp.real.matrix<-function(A,eigenSigns=NULL,mCSigns=NULL,tol=1e-15){
    n<-nrow(A)    
    vParams<-rep(NA,n^2)
    if (ncol(A)!=n){print("A must be a square matrix")}
    else{
	if (n==1){vParams<-A[1,1];mCSigns=1}
	else{
	    eigA<-eigen(A)
	    if (!is.null(eigenSigns)){
		eigA$vectors<-apply(rbind(eigenSigns,eigA$vectors),2,function(x){if (sign(x[2])!=x[1]){x<-(-1)*x};x[-1]}) ## we want to fix the eigenvector matrix
	    }
	    if(is.complex(eigA$values)){
		print("A has complex eigen values, calling complex parametrization!")
		vParams<-.par.inv.transform.decomp.matrix(A,eigenSigns,mCSigns)
	    }else{
	        k<-n*(n-1)/2
		vParams[1:n]<-eigA$values
		qrP<-qr(eigA$vectors)   		
		QqrP<-qr.Q(qrP)
		RqrP<-qr.R(qrP)
		vNegR<-which(diag(RqrP<0))
		RqrP[vNegR,]<-RqrP[vNegR,]*(-1)
		QqrP[,vNegR]<-QqrP[,vNegR]*(-1)
		lQparams<-.par.inv.transform.orth.matrix.givens(QqrP,mCSigns,tol)
		vParams[(n+1):(n+k)]<-lQparams$vParams
		mCSigns<-lQparams$cSigns
		RqrP<-sapply(1:ncol(RqrP),function(i,mR){
		    k<-nrow(mR)
		    orgcolR<-mR[,i]
		    colR<-c(rep(NA,i-1),rep(0,n-i+1))
		    if (i>1){## This turned out to be slightly faster then a while loop
			if (i>2){colR[2:(i-1)]<-tan(orgcolR[2:(i-1)]*0.5*pi/(sqrt(1-cumsum(orgcolR[1:(i-2)]^2))))}
			colR[1]<-tan(orgcolR[1]*0.5*pi)		    
		    }
		    colR
		},mR=RqrP,simplify=TRUE)
		vParams[(n+k+1):(n^2)]<-RqrP[upper.tri(RqrP)]
	    }	    
	}
    }
    list(vParams=vParams,cSigns=mCSigns)
}
## -------------------------------------------------------------------------------------------------
## ----------------- general case ------------------------------------------------------------------
.par.transform.decomp.matrix<-function(vParams,n){
    print("Cannot do parametrization with general complex eigenvalues yet")
    matrix(NA,n,n)
}

.par.inv.transform.decomp.matrix<-function(A,eigenSigns,mCSigns){
    print("Cannot do deparametrization with general complex eigenvalues yet")
    rep(NA,nrow(A)^2)
}
## -------------------------------------------------------------------------------------------------
## -------------------------------------------------------------------------------------------------


## ----------- parametrization of orthogonal matrix ------------------------------------------------
.par.transform.orth.matrix.givens<-function(vParams,n,mCSigns=NULL,vDiagQ=NULL){
    Q<-NA
    k<-n*(n-1)/2
    if (is.null(vDiagQ)){vDiagQ<-rep(1,n)}
    if (is.null(mCSigns)){mCSignsQ<-matrix(1,2,k)}
    else{mCSignsQ<-mCSigns}
    if (length(vParams)!=k){print(paste("Cannot calculate Q, vParams should be of lenght : ",n*(n-1)/2));}
    else{
	if (n==1){Q<-matrix(vParams[1],1,1)}## vParams[1] has to be 1 or -1 for orthogonality
	else{
	    Q<-diag(vDiagQ,n,n)
	    ij<-k
    	    for (j in (n-1):1){
		for (i in n:(j+1)){
		    rho<-vParams[ij]
		    if (abs(rho)<1.5){rho<-rho/3}
		    else{rho<-rho*1.5}
		    cs<-.inv.rho.givens.real(rho)
		    if(!is.null(mCSigns)){
			if (sign(cs[1])!=mCSigns[1,ij]){cs[1]<-(-1)*cs[1] }
			if (sign(cs[2])!=mCSigns[2,ij]){cs[2]<-(-1)*cs[2] }
		    }
		    else{mCSignsQ[,ij]<-sign(cs)}
		    Q<-.givens.mult.M(Q,cs[1],-cs[2],j,i) 
		    ij<-ij-1;
		}
	    }
	}
    }
    list(Q=Q,cSigns=mCSignsQ,vDiagQ=vDiagQ)
}

.par.inv.transform.orth.matrix.givens<-function(Q,mCSigns=NULL,tol=1e-15,check1=FALSE){
    n<-nrow(Q)
    k<-n*(n-1)/2
    vParams<-rep(NA,k)
    if (is.null(mCSigns)){mCSignsQ<-matrix(1,2,k)}
    else{mCSignsQ<-mCSigns}
    if (ncol(Q)!=n){print("Q must be a square matrix")}
    else{
	if (n==1){vParams<-Q[1,1]}
	else{
	    ij<-1
	    GtQ<-Q
	    for (j in 1:(n-1)){
		for (i in (j+1):n){
		    cs<-.givens.real(GtQ[j,j],GtQ[i,j],getrho=FALSE)
		    if (!is.null(mCSigns)){
			if (sign(cs[1])!=mCSigns[1,ij]){cs[1]<-(-1)*cs[1]}
			if (sign(cs[2])!=mCSigns[2,ij]){cs[2]<-(-1)*cs[2]}		
		    }else{mCSignsQ[,ij]<-sign(cs)}
		    GtQ<-.givens.mult.M(GtQ,cs[1],cs[2],j,i) ## j<i so it should be given in this order so that rows not swapped
		    rho<-.givens.real(NULL,NULL,TRUE,cs)				    
		    if (abs(rho)<0.5){rho<-rho*3}
		    else{rho<-rho/1.5}
		    vParams[ij]<-rho
		    ij<-ij+1
		}
	    }
	    vDiagGtQ<-diag(GtQ)
	    ## this is a check for whether Q is orthogonal if it is we should have an identity matrix in GtQ (up to numerical prescision)
	    if (check1){vParams<-c(sign(diag(GtQ))*rep(1,n),vParams)}## add the sign of the diagonal	    
	    if (!is.null(tol)){## check whether we have calculated ok
	        diag(GtQ)<-abs(diag(GtQ)) ## we might have +/- 1 on the diagonal 
		GtQ[which(abs(GtQ)<tol)]<-0
		vCorrectNum<-which(abs(abs(GtQ)-diag(1,n,n))<tol) ## if on the diagonal
		GtQ[vCorrectNum]<-diag(1,n,n)[vCorrectNum]
		invGtQ<-solve(GtQ)
		if (length(GtQ[!GtQ==invGtQ])>0){print(paste("Q matrix was not orthogonal or perhaps numerical errors, errors : ",GtQ[!GtQ==invGtQ]-invGtQ[!GtQ==invGtQ]))}
	    }
	}
    }
    list(vParams=vParams,cSigns=mCSignsQ,vDiagQ=vDiagGtQ)
}
## -------------------------------------------------------------------------------------------------

## ----------------- Transform vector to positive numbers ------------------------------------------
.num.transform.positive<-function(x,method){
    y<-0
    if (method=="exp"){y<-exp(x)}
    if (method=="atan"){y<-atan(x)+pi/2}
    y
}

.num.inv.transform.positive<-function(y,method){
    x<-0
    if (method=="exp"){x<-log(y)}
    if (method=="atan"){x<-tan(y-pi/2)}
    x
}

## -------------------------------------------------------------------------------------------------

## ----------- Some additional transformations -------------------------------------------------
#.par.transform.normal.matrix<-function(vParams,kY){
#    eigenM<-diag(vParams[1:kY],ncol=kY,nrow=kY)
#    P<-.par.transform.orth.matrix(vParams[(kY+1):length(vParams)])
#    P%*%eigenM%*%solve(P)
#}

#.par.inv.transform.normal.matrix<-function(A){
#    eigA<-eigen(A)
#    c(eigA$values,.par.inv.transform.orth.matrix(eigA$vectors))
#}

.par.transform.orth.matrix.cayley<-function(vParams){
## This is Cayley's parametrization of an orthognal matrix as found on PlanetMath.org
## vParams represent a skew-symmetric matrix
## The only restriction on the orthogonal matrix P is that it does not have
## an eigenvalue equal to -1, which should in our case be OK
## as we're estimating so there is a prob of 0 of getting exact value
    kY<-(1+sqrt(1+8*length(vParams)))/2
    CayleyA<-matrix(0,nrow=kY,ncol=kY)
    CayleyA[upper.tri(CayleyA)]<-vParams
    CayleyA[lower.tri(CayleyA)]<-(-1)*vParams
    (diag(1,nrow=kY,ncol=kY)+CayleyA)%*%solve(diag(1,nrow=kY,ncol=kY)+CayleyA)
}

.par.inv.transform.orth.matrix.cayley<-function(A){
## This is Cayley's parametrization of an orthognal matrix as found on PlanetMath.org
## vParams represent a skew-symmetric matrix
## The only restriction on the orthogonal matrix P is that it does not have
## an eigenvalue equal to -1, which should in our case be OK
## as we're estimating so there is a prob of 0 of getting exact value
    CayleyA<-solve(A+diag(1,nrow=nrow(A),ncol=ncol(A)))%*%(A-diag(1,nrow=nrow(A),ncol=ncol(A)))
    CayleyA[upper.tri(CayleyA)] ## this is done bycol, upper.tri does not have option byrow=TRUE ...
}
## -------------------------------------------------------------------------------------------------
