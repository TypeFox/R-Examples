# KR <-
# function(object, ...) UseMethod('KR')

KR.varComp=function(object, Lmat, Vbet, svd.VLbet, X, K, V, ...)
{
	bhat=coef(object, 'fixed'); 
	if(length(bhat)== 0L){
		return(structure(numeric(0L), names=character(0L), numDF = numeric(0L), Scale=numeric(0L), `F value`=numeric(0L), `Pr(>F)`=numeric(0L), vcov.beta = matrix(NA_real_, 0L, 0L)))
	}
	#if(missing(K)) 
		K= model.matrix(object, what='K')
	if(missing(X)) X= model.matrix(object, what='X')
  if(missing(Lmat)) {Lmat=diag(1, ncol(X)); rownames(Lmat)=colnames(X)}
  if(!is.matrix(Lmat))	Lmat = matrix(Lmat, 1L, dimnames = list("", names(Lmat)))
  
  if(!is.null(colnames(Lmat)) && !is.null(colnames(X))){
	if(ncol(Lmat) < ncol(X)){
		L=matrix(0, nrow(Lmat), ncol(X))
		colnames(L) = colnames(X)
		rownames(L) = rownames(Lmat)
		L[, colnames(Lmat)]=Lmat
		Lmat=L
	}else if (ncol(Lmat) == ncol(X)){
		Lmat = Lmat[, colnames(X), drop=FALSE]
	}else stop("`Lmat` has more columns than the number of fixed effect parameters.")	
  }
  if(ncol(Lmat) != ncol(X)) stop("`Lmat` has incorrect number of columns.")
  if(is.null(rownames(Lmat))) rownames(Lmat) = rep('', nrow(Lmat))
  if(missing(V)) V= vcov(object, 'Y'); VI=solve(V)
	if(missing(Vbet)) Vbet=vcov(object, what='beta', beta.correction=FALSE)
	Phi = Vbet
	p=length(bhat)
	ell=qr(Lmat)$rank
	r=length(coef(object, 'varComp'))
	J=r-1
	n=nrow(V)
	
	w=vcov.varComp(object, 'varComp', drop=FALSE)
	w=as.matrix(nearPD(w)$mat)
	tmp=which(object$parms==0)
	w[tmp,]=w[,tmp]=0

	
	VIX=solve(V, X)
	# Phi = solve(crossprod(X, VIX))

	Q=array(NA_real_, dim=c(p,p,r,r))
		for(i in seq_len(J)){
			for(j in seq_len(J))
				Q[,,i,j]=crossprod(VIX, K[[i]]%*%VI%*%K[[j]]%*%VIX)
			Q[,,i,r]=crossprod(VIX, K[[i]]%*%VI%*%VIX)
		}
		for(j in seq_len(J)) Q[,,r,j]=crossprod(VIX, VI%*%K[[j]]%*%VIX)
		Q[,,r,r]=crossprod(VIX, VI%*%VIX)
	P=array(NA_real_, dim=c(p,p,r))
		for(i in seq_len(J)) P[,,i]=-crossprod(VIX, K[[i]]%*%VIX)
		P[,,r]=-crossprod(VIX)

	LambdaTilde=matrix(0, p,p)
		for(i in seq_len(r))
			for(j in seq_len(r))
				LambdaTilde=LambdaTilde + w[i,j]*(Q[,,i,j]-P[,,i]%*%Phi%*%P[,,j])
		LambdaTilde=Phi%*%LambdaTilde%*%Phi

	Rstar=0
	Phi_A = Vbet + 2 * LambdaTilde - Rstar

	if(ell==0L) return(structure(NA_real_, names='denDF', numDF = ell, Scale=NA_real_, `F value`=NA_real_, `Pr(>F)`=NA_real_, vcov.beta = Phi_A))
	
	Lmat=t(Lmat)
	F=drop( crossprod(bhat, Lmat%*%solve(crossprod(Lmat, Phi_A%*%Lmat))%*%crossprod(Lmat, bhat))) / ell
	Theta = Lmat%*%solve(crossprod(Lmat, Phi%*%Lmat), t(Lmat))
	PTP=Phi%*%Theta%*%Phi
	
	A1=0;	for(i in seq_len(r)) for(j in seq_len(r)) A1=A1+w[i,j]*sum(PTP*P[,,i])*sum(PTP*P[,,j])
	A2=0;	for(i in seq_len(r)) for(j in seq_len(r)) A2=A2+w[i,j]*sum(diag(PTP%*%P[,,i]%*%PTP%*%P[,,j]))
	B=(A1+6*A2)/2/ell
	
	g=((ell+1)*A1-(ell+4)*A2)/(ell+2)/A2; c=3*ell+2*(1-g); 
	d1=g/c; d2=(ell-g)/c; d3=(ell-g+2)/c
	Estar=1/(1-A2/ell); Vstar=2/ell*(1+d1*B)/(1-d2*B)^2/(1-d3*B)
	rhoTilde=Vstar/2/Estar/Estar
	m=4+(2+ell)/(ell*rhoTilde-1); lambda=m/Estar/(m-2)
	
	structure(m, names='denDF', numDF = ell, Scale=lambda, `F value`=F*lambda, `Pr(>F)`=pf(lambda*F, ell, m, lower.tail=FALSE), vcov.beta = Phi_A)
}

