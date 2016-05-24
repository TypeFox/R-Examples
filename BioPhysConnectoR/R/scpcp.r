#############################################
#   This code is subject to the license as stated in DESCRIPTION.
#   Using this software implies acceptance of the license terms:
#    - GPL 2
#
#   (C) by F. Hoffgaard, P. Weil, and K. Hamacher in 2009.
#
#  keul(AT)bio.tu-darmstadt.de
#
#
#  http://www.kay-hamacher.de
#############################################

scpcp<-function(T,R,cm,pstart=0.5,maxiter=2000,chains=NULL,maxtol=1e-11,file=NULL,im=NULL){

	compute_B_and_G<-function(n, B){
		prefac<-8*pi*pi/(n*n)
		ms1<-n+1
		b<-rep(0,n)
		g<-matrix(0,n,n)
		if(!is.null(file)){
		  cat("G-Matrix:\n",file=file,append=TRUE)
		}
		for(i in 1:n){
			b[i]<-(n-1)*(n-1)*B[i,i]
			z<-0
			for(k in 1:n){
				if(i!=k){
					z<-z+B[k,k]-2*(n-1)*B[i,k]
					for(ks in 1:n){
						if(i!=ks && k!=ks){
							z<-z+B[k,ks]
						}
					}
				}
			}
			b[i]<-b[i]+z
			b[i]<-b[i]*prefac
		}
		
		for(i in 1:n){
			for(k in 1:n){
				g[i,k]<-3*(B[i,i]+B[k,k]-B[i,k]-B[k,i])
				if(!is.null(file)){
					cat(i,k,g[i,k],"\n",file=file,append=TRUE)
				}
			}
		}
		return(list(b=b,g=g))
	}

	compute_Pi<-function(n,rc,cm,p,covind){
		diag(cm)<-diag(cm[,-1])<-diag(cm[-1,])<-0
		p<-cm*p
		C<-colSums(cm)
		z<-colSums(p)
		Pi<-z/C
		Pi[which(C==0)]<-1
		return(Pi)
	}

	compute_Q<-function(rc, cm, p){
		z<-0
		C<-0
		cm2<-cm
		diag(cm2)<-diag(cm2[,-1])<-diag(cm2[-1,])<-0
		C2<-sum(cm2)
		z2<-sum(cm2*p)
		return(z2/C2)
	}

	compute_F<-function(T, R, n, rc, cm, p, im, m){
		F<- -1.5*T*n*log(2*pi)
		cm2<-cm
		diag(cm2[,-1])<-diag(cm2[-1,])<-0
		z2<-sum(cm2*p*im)
		F<-F - R*R*z2/2
		E<-3*(n-1)*T/2 - R*R*z2/2
		det<-compute_detL(m=m)
		F<-F - 1.5*T*log(abs(det))
		F<-F + 1.5*T*log(2*pi/10)
		return(list(F=F,E=E,det=det))
	}

	compute_detL<-function(m){
		oo<-get.svd(m)$ev
		detU<-1/prod(oo)
		return(detU)
	}

	initM<-function(rc, T2, beta2, cm, p, n, im, XI){
		M2<-matrix(0,n,n)
		M<-matrix(0,n,n)
		rc2<-rc-1	
		XI2<-matrix(0,n,n)
		for(i in 1:n){
			for(j in 1:n){
				if(XI[i,j]==TRUE){
					XI2[i,j]<-1
				}
			}
		}
		nrc<-length(rc[,1])
		out<-.C("initM",r=as.integer(rc2[,1]),c=as.integer(rc2[,2]),n=as.integer(n),nrc=as.integer(nrc),XI=as.integer(XI),im=as.double(im),p=as.double(p),T2=as.double(T2),beta2=as.double(beta2),M=as.double(M2),PACKAGE="BioPhysConnectoR")
		M2<-matrix(out$M,ncol=n)
		return(M2)
	}

	init_p_with_cm<-function(value,rc,cm){
		p2<-matrix(0,dim(cm)[1],dim(cm)[2])
		p2<-cm*value
		diag(p2)<-0
		return(p2)
	}

	contact_in_p<-function(p,covind){
		diag(p[(-1),])[covind]<-1
		diag(p[,(-1)])[covind]<-1
		return(p)
	}

#	prep of in same-chain matrix XI and index for all covalently
#	bound contacts (off-diagonal)
	dcm<-dim(cm)
	n<-dcm[1]
	XI<-matrix(0,dcm[1],dcm[2])
	if(is.null(chains)){
		covind<-1:(n-1)
		XI<-matrix(1,dcm[1],dcm[2])
	}else{
		chainEnd<-cumsum(chains)
		chainStart<-c(1,chainEnd+1)
		chainStart<-chainStart[-length(chainStart)]
		nchains<-length(chainEnd)
		covind<-1:(n-1)
		covind<-covind[-chainEnd[-nchains]]
		for(i in 1:nchains){
			a<-chainStart[i]
			b<-chainEnd[i]
			XI[a:b,a:b]<-1
		}
	}
#	prep of contact matrix to list
	rc<-which(cm==1,arr.ind=TRUE)
	nrc<-length(rc[,1])
#	cont preparation
	R22<-(R^2)/2
	T2<-T/2
	beta2<-2/T
#	initiate
	p<-init_p_with_cm(value=pstart,rc=rc,cm=cm)
	p<-contact_in_p(p=p,covind=covind)

#	iterate
	for(iter in 1:maxiter){
		tol<-0
		m<-initM(rc=rc,T2=T2,beta2=beta2,cm=cm,p=p,n=n,im=im,XI=XI)
		B<-matrix.inverse(m)
		f<-compute_F(T=T,R=R,n=n,rc=rc,cm=cm,p=p,im=im,m=m)
		e<-f$E
		det<-f$det
		f<-f$F
		error<-0
		for(i in 1:dim(rc)[1]){
			ri<-rc[i,1]
			ci<-rc[i,2]
			if(ri!=ci && !(XI[ri,ci]==1 && (ci==ri+1 || ri==ci+1)) ){
				Gij<-(B[ri,ri]+B[ci,ci]-B[ci,ri]-B[ri,ci])
				if(Gij<=0){
					pij<-1
				}else{
					pij<-(1-pgamma(R22/Gij,shape=3/2,lower.tail=FALSE))
				}
				error<-abs(pij-p[ri,ci])
				if(error>tol){
					tol<-error
				}
				p[ri,ci]<-pij
			}
		}
		p<-contact_in_p(p=p,covind=covind)
		q<-compute_Q(rc=rc,cm=cm,p=p)
		if(!is.null(file)){
			cat("\nf: ",f,"\te: ",e,"\tQ: ",q," \ttol: ",round(tol,digits=8),"\terr: ",error,"\n",sep="",file=file,append=TRUE)
		}else{
			cat("\nf: ",f,"\te: ",e,"\tQ: ",q," \ttol: ",round(tol,digits=8),"\terr: ",error,"\n",sep="")
		}
		if(tol<maxtol){
			break
		}
	}
	bg<-compute_B_and_G(n=n,B=B)
	bfacs<-bg$b
	gmat<-bg$g
	entropy<-(e-f)/T
	Pi<-compute_Pi(n=n,rc=rc,cm=cm,p=p,covind=covind)
	if(!is.null(file)){
		cat("Bfacs:",bfacs,sep="\n",file=file,append=TRUE)
		cat("Pi:",Pi,sep="\n",file=file,append=TRUE)
	}
	return(list(free=f,intern=e,entropy=entropy,q=q,bfacs=bfacs,pi=Pi,gmat=gmat,iter=iter,err=error))
}
