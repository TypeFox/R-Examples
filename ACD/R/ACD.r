#ACD: categorical data analysis with complete or missing responses
#Date: 2012/11/23
#Authors: Frederico Zanqueta Poleto, Julio da Motta Singer, Carlos Daniel Paulino, Fabio Mathias Correa and 
#  Enio Galinkin Jelihovschi
#License: The functions may be distributed free of charge and used by anyone if credit is given. It has been tested
#         fairly well, but it comes with no guarantees and the authors assume no liability for its use or misuse.
#URL: http://www.poleto.com/missing.html



#This function receives the input of categorical data for analysis of complete or missing categorical data,
#checks some consistencies and creates some of the quantities that the other functions will use in the subsequent analyses.

readCatdata<-function(TF=NULL, Zp=NULL, Rp=NULL) {
	if(is.null(TF)) {
		stop("Specify the table of frequencies (matrix TF).")
	}
	if(!all(TF>=0)) {
		stop("The table of frequencies (matrix TF) must have only non-negative numbers.")
	}
	if(!is.matrix(TF)) {
		TF<-t(TF)
	}
	S<-nrow(TF)
	Tt<-rep(1,S)

	if( (is.null(Zp)) && (is.null(Rp)) ) {

		miss<-FALSE
		R<-ncol(TF)

		Rp<-matrix(0,S,1)
		l<-rep(0,S)

		Zs<-vector("list",S) #\{Z_s\}
		Zbs<-vector("list",S) #\{\bar{Z}_s\}
		Ns<-vector("list",S) #\{N_s\}
		Nsm<-vector("list",S) #\{N_{s+}\}
		ps<-vector("list",S) #\{p_s\}
		pbs<-vector("list",S) #\{\bar{p}_s\}
		Vps<-vector("list",S) #\{\Sigma(p_s)\}
		Vpbs<-vector("list",S) #\{\Sigma(\bar{p}_s)\}
		Sn<-vector("character",S)

		Nst<-vector("list",sum(Tt)) #\{N_s\}
		Nbst<-vector("list",sum(Tt)) #\{N_s\}
		pst<-vector("list",sum(Tt)) #\{p_s\}
		pbst<-vector("list",sum(Tt)) #\{\bar{p}_s\}
		Vpst<-vector("list",sum(Tt)) #\{\Sigma(p_{st})\}
		Vpbst<-vector("list",sum(Tt)) #\{\Sigma(\bar{p}_{st})\}
		STSn<-vector("character",sum(Tt))

		nstm<-matrix(0,S,max(Tt)) #\{n_{st+}\}

		nsmm<-rowSums(TF) #\{n_{s++}\}
		nstm<-as.matrix(nsmm)
		theta<-numeric(S*R)
		Vtheta<-matrix(0,S*R,S*R)
		for(s in 1:S) {
			Zs[[s]]<-cbind(diag(R))
			Zbs[[s]]<-cbind(diag(R-1))
			Nsm[[s]]<-rep(nsmm[s],R)
			Nst[[s]]<-Ns[[s]]<-TF[s,]
			Nbst[[s]]<-TF[s,1:(R-1)]
			pst[[s]]<-ps[[s]]<-theta[((s-1)*R+1):(s*R)]<-TF[s,]/nsmm[s]
			pbst[[s]]<-pbs[[s]]<-ps[[s]][1:(R-1)]
			Vpst[[s]]<-Vps[[s]]<-Vtheta[((s-1)*R+1):(s*R),((s-1)*R+1):(s*R)]<-(diag(ps[[s]])-ps[[s]]%*%t(ps[[s]]))/nsmm[s]
			Vpbst[[s]]<-Vpbs[[s]]<-Vps[[s]][1:(R-1),1:(R-1)]
			Sn[s]<-as.character(paste("s",s,sep=""))
			STSn[s]<-as.character(paste("st",s,".",1,sep=""))
		}
		names(Zs)<-names(Zbs)<-names(Ns)<-names(Nsm)<-names(ps)<-names(pbs)<-names(Vps)<-names(Vpbs)<-Sn
		names(Nst)<-names(Nbst)<-names(pst)<-names(pbst)<-names(Vpst)<-names(Vpbst)<-STSn

		thetab<-c( kronecker(diag(S),cbind(diag(R-1),rep(0,R-1)))%*%theta )
		Vthetab<-kronecker(diag(S),cbind(diag(R-1),rep(0,R-1)))%*%Vtheta%*%t(kronecker(diag(S),cbind(diag(R-1),rep(0,R-1))))
		bs<-c(rep(0,R-1),1) #b_s
		Bs<-rbind(diag(R-1),rep(-1,R-1)) #B_s
		b<-kronecker(rep(1,S),bs) #b
		B<-kronecker(diag(S),Bs) #B

		res<-list(TF=TF,R=R,S=S,call=match.call(),miss=miss,Rp=Rp,l=l,Tt=Tt,Zs=Zs,Zbs=Zbs,
			  Ns=Ns,Nst=Nst,Nbst=Nbst,Nsm=Nsm,nstm=nstm,nsmm=nsmm,ps=ps,pst=pst,pbs=pbs,pbst=pbst,Vps=Vps,Vpbs=Vpbs,Vpst=Vpst,Vpbst=Vpbst,
			  bs=bs,Bs=Bs,b=b,B=B,theta=theta,Vtheta=Vtheta,thetab=thetab,Vthetab=Vthetab)
	} else {

		if(is.null(Rp) | is.null(Zp)) {
			stop("Specify Zp AND Rp if you have missing data or remove both matrices if you have complete data.")
		}
		miss<-TRUE
		R<-nrow(Zp)
		if(!is.matrix(Rp)) {
			Rp<-t(Rp)
		}
		if(!all((Rp>=2) | (Rp==0))) {
			stop("The matrix Rp must have only numbers >=2 or equal 0.")
		}
		if(!all((Zp==0) | (Zp==1))) {
			stop("The values in the matrix Zp must be equal 1 or 0.")
		}
		l<-rowSums(Rp) #\{l_s\}
		if(ncol(Zp)!=sum(l)) {
			stop(paste("The number of columns of Zp should be equal to ",sum(l)," according to your Rp.",sep=""))
		}
		if(ncol(TF)!=(R+max(l))) {
			stop(paste("The number of columns of TF should be equal to ",R+max(l)," according to your Zp and Rp.",sep=""))
		}
		for(s in 1:S) {
			for(i in 1:ncol(Rp)) {
				if((Tt[s]==i) && (Rp[s,i]>0)) {
					Tt[s]<-Tt[s]+1
				}
			}
		}

		Zs<-vector("list",S) #\{Z_s\}
		Zbs<-vector("list",S) #\{\bar{Z}_s\}
		Ns<-vector("list",S) #\{N_s\}
		Nsm<-vector("list",S) #\{N_{s+}\}
		ps<-vector("list",S) #\{p_s\}
		pbs<-vector("list",S) #\{\bar{p}_s\}
		Sn<-vector("character",S)

		Zst<-vector("list",sum(Tt)-S) #\{Z_{st}\}
		Zbst<-vector("list",sum(Tt)-S) #\{\bar{Z}_{st}\}
		STn<-vector("character",sum(Tt)-S)

		Nst<-vector("list",sum(Tt)) #\{N_s\}
		Nbst<-vector("list",sum(Tt)) #\{N_s\}
		pst<-vector("list",sum(Tt)) #\{p_s\}
		pbst<-vector("list",sum(Tt)) #\{\bar{p}_s\}
		Vpst<-vector("list",sum(Tt)) #\{\Sigma(p_{st})\}
		Vpbst<-vector("list",sum(Tt)) #\{\Sigma(\bar{p}_{st})\}
		STSn<-vector("character",sum(Tt))

		nstm<-matrix(0,S,max(Tt)) #\{n_{st+}\}

		i1<-0
		i2<-1
		i3<-0
		for(s in 1:S) {
			if(Tt[s]>1) {
				Zs[[s]]<-cbind(diag(R),Zp[,i2:(i2+sum(Rp[s,1:(Tt[s]-1)])-1)])
			} else {
				Zs[[s]]<-diag(R)
			}
			Ns[[s]]<-TF[s,1:(R+l[s])]
			Sn[s]<-as.character(paste("s",s,sep=""))
			i3<-i3+1
			Nst[[i3]]<-TF[s,1:R]
			Nbst[[i3]]<-TF[s,1:(R-1)]
			Nsm[[s]]<-rep(sum(Nst[[i3]]),R)
			nstm[s,1]<-sum(Nst[[i3]])
			ps[[s]]<-pst[[i3]]<-Nst[[i3]]/nstm[s,1]
			pbs[[s]]<-pbst[[i3]]<-Nbst[[i3]]/nstm[s,1]
			Vpst[[i3]]<-(diag(pst[[i3]])-pst[[i3]]%*%t(pst[[i3]]))/nstm[s,1]
			Vpbst[[i3]]<-Vpst[[i3]][1:(R-1),1:(R-1)]
			STSn[i3]<-as.character(paste("st",s,".1",sep=""))
			zz<-diag(R-1)
			i4<-R+1
			if(Tt[s]>1) {
				for(tt in 2:Tt[s]) {
					i1<-i1+1
					i3<-i3+1
					Zst[[i1]]<-Zp[,i2:(i2+Rp[s,(tt-1)]-1)]
					if(!all(  rowSums( Zst[[i1]] )==rep(1,R)  )) {
						stop(paste("The submatrix of Zp associated with s=",s,",t=",tt," is not a partition of the hypothetic completed categorized data, because the sum of the columns vectors is not equal a vector of 1's.",sep=""))
					}
					Zbst[[i1]]<-Zp[1:(R-1),i2:(i2+Rp[s,(tt-1)]-2)]
					zz<-cbind(zz,Zbst[[i1]])
					STSn[i3]<-STn[i1]<-as.character(paste("st",s,".",tt,sep=""))
					i2<-i2+Rp[s,(tt-1)]
					Nst[[i3]]<-TF[s,i4:(i4+Rp[s,(tt-1)]-1)]
					Nbst[[i3]]<-TF[s,i4:(i4+Rp[s,(tt-1)]-2)]
					Nsm[[s]]<-c(Nsm[[s]],rep(sum(Nst[[i3]]),Rp[s,(tt-1)]))
					nstm[s,tt]<-sum(Nst[[i3]])
					pst[[i3]]<-Nst[[i3]]/nstm[s,tt]
					ps[[s]]<-c(ps[[s]],pst[[i3]])
					pbst[[i3]]<-Nbst[[i3]]/nstm[s,tt]
					pbs[[s]]<-c(pbs[[s]],pbst[[i3]])
					Vpst[[i3]]<-(diag(pst[[i3]])-pst[[i3]]%*%t(pst[[i3]]))/nstm[s,tt]
					Vpbst[[i3]]<-Vpst[[i3]][1:(Rp[s,(tt-1)]-1),1:(Rp[s,(tt-1)]-1)]
					i4<-i4+Rp[s,(tt-1)]
				}
			}
			Zbs[[s]]<-zz
		}

		names(Zs)<-names(Zbs)<-names(Ns)<-names(Nsm)<-names(ps)<-names(pbs)<-Sn
		names(Zst)<-names(Zbst)<-STn
		names(Nst)<-names(Nbst)<-names(pst)<-names(pbst)<-names(Vpst)<-names(Vpbst)<-STSn

		nsmm<-rowSums(TF) #\{n_{s++}\}
		bs<-c(rep(0,R-1),1) #b_s
		Bs<-rbind(diag(R-1),rep(-1,R-1)) #B_s
		b<-kronecker(rep(1,S),bs) #b
		B<-kronecker(diag(S),Bs) #B

		res<-list(TF=TF,R=R,S=S,call=match.call(),miss=miss,Rp=Rp,l=l,Tt=Tt,Zs=Zs,Zst=Zst,Zbs=Zbs,Zbst=Zbst,
			  Ns=Ns,Nst=Nst,Nbst=Nbst,Nsm=Nsm,nstm=nstm,nsmm=nsmm,ps=ps,pst=pst,pbs=pbs,pbst=pbst,Vpst=Vpst,Vpbst=Vpbst,
			  bs=bs,Bs=Bs,b=b,B=B)

	}
	class(res)<-"readCatdata"
	invisible(res)
}



print.readCatdata<-function(x, digits=max(3,getOption("digits")-3), ...) {
	cat("\nCall: ",deparse(x$call),"\n",sep="")

	cat("\nS=",x$S," subpopulations  x  R=",x$R," response categories  with ",ifelse(x$miss,"MISSING","COMPLETE")," data\n",sep="")

	pp<-matrix(0,x$S,x$R)
	epp<-matrix(0,x$S,x$R)
	for(s in 1:x$S) {
		if(x$miss) {
			pp[s,]<-x$pst[[paste("st",s,".",1,sep="")]]
			epp[s,]<-sqrt(diag(x$Vpst[[paste("st",s,".",1,sep="")]]))
		} else {
			pp[s,]<-x$ps[[paste("s",s,sep="")]]
			epp[s,]<-sqrt(diag(x$Vps[[paste("s",s,sep="")]]))
		}
	}

	cat("\nProportions",ifelse(x$miss," of the complete data",""),":\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	cat("\nStandard errors of the proportions",ifelse(x$miss," of the complete data",""),":\n",sep="")
	print.default(format(round(epp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	cat("\n")
}



summary.readCatdata<-function(object, digits=max(3,getOption("digits")-3), ...) {

	cat("\nCall: ",deparse(object$call),"\n\n",sep="")
	cat("\nS=",object$S," subpopulations  x  R=",object$R," response categories  with ",ifelse(object$miss,"MISSING","COMPLETE")," data\n",sep="")

	cat("\n\nTable of frequencies",ifelse(object$miss," of the complete data",""),":\n",sep="")
	print.default(object$TF[1:object$S,1:object$R],quote=FALSE,print.gap=2)

	pp<-matrix(0,object$S,object$R)
	epp<-matrix(0,object$S,object$R)
	for(s in 1:object$S) {
		if(object$miss) {
			pp[s,]<-object$pst[[paste("st",s,".",1,sep="")]]
			epp[s,]<-sqrt(diag(object$Vpst[[paste("st",s,".",1,sep="")]]))
		} else {
			pp[s,]<-object$ps[[paste("s",s,sep="")]]
			epp[s,]<-sqrt(diag(object$Vps[[paste("s",s,sep="")]]))
		}
	}

	cat("\nProportions",ifelse(object$miss," of the complete data",""),":\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	cat("\nStandard errors of the proportions",ifelse(object$miss," of the complete data",""),":\n",sep="")
	print.default(format(round(epp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	if(object$miss) {
		cat("\n\nMissing data frequencies and associated column vectors indicating\nthe relation with the original set of R response categories:\n",sep="")
		for(s in 1:object$S) {
			if(object$Tt[s]>1) {
				if(object$S>1) {
					cat("\nSubpopulation ",s,sep="")
				}
				for(tt in 2:object$Tt[s]) {
					cat("\n")
					zz<-t(rbind(object$Nst[[paste("st",s,".",tt,sep="")]], rep("",object$Rp[s,tt-1]), rep("",object$Rp[s,tt-1]), object$Zst[[paste("st",s,".",tt,sep="")]]))
					colnames(zz)<-c(rep("",3),paste("[,",1:object$R,"]",sep=""))
					print.default(zz,quote=FALSE,print.gap=2,right=TRUE)
				}
			}
		}
	}
	cat("\n")
}





#Fit saturated models for the marginal probabilities of categorization under MAR/MCAR by maximum likelihood
satMarML<-function(catdataobj,missing="MAR",method="EM",start=NULL,zero=NULL,maxit=100,trace=0,epsilon1=1e-6,epsilon2=1e-6,zeroN=NULL, digits=max(3,getOption("digits")-3)) {
	if(class(catdataobj)!="readCatdata") {
		stop("catdataobj is not of the class readCatdata")
	}
	tole<-1e-323 #tol argument for solve function

	R<-catdataobj$R
	S<-catdataobj$S
	miss<-catdataobj$miss
	Rp<-catdataobj$Rp
	l<-catdataobj$l
	Tt<-catdataobj$Tt
	Zs<-catdataobj$Zs
	Zst<-catdataobj$Zst
	Zbs<-catdataobj$Zbs
	Zbst<-catdataobj$Zbst
	Ns<-catdataobj$Ns
	Nst<-catdataobj$Nst
	Nsm<-catdataobj$Nsm
	nstm<-catdataobj$nstm
	nsmm<-catdataobj$nsmm
	ps<-catdataobj$ps
	pst<-catdataobj$pst
	pbst<-catdataobj$pbst
	bs<-catdataobj$bs
	Bs<-catdataobj$Bs
	b<-catdataobj$b
	B<-catdataobj$B

	if(miss==FALSE) {
		stop("This function is just to categorical data with missings.\ntheta in your object already contains the ML estimates of the probabilities.")
	}
	if(!any(missing==c("MAR","MCAR"))) {
		stop("missing (mechanism) can be only \"MAR\" or \"MCAR\" ")
	}
	if(!any(method==c("EM","FS-MCAR","NR/FS-MAR"))) {
		stop("method (of ML estimation) can be only \"EM\", \"FS-MCAR\" or \"NR/FS-MAR\" ")
	}
	if(!is.null(zero)) {
		if(any(zero<0 | zero>0.5)) {
			stop("The elements of the vector zero must be nonnegative and less than or equal 0.5.")
		}
		if((S>1) && (length(zero)<S)) {
			zero<-c(zero,rep(zero[1],S-length(zero)))
		}
	} else {
		zero<-1/(R*nstm[,1])
	}
	if(!is.null(zeroN)) {
		if(!is.matrix(zeroN) || (nrow(zeroN)!=S) || (ncol(zeroN)!=max(Tt))) {
			if(is.vector(zeroN,mode="numeric")) {
				zeroN<-matrix(zeroN[1],S,max(Tt))
			} else {
				stop("zeroN must be a number or a matrix with dimensions S x (max(Tt))")
			}
		}
		if(any(zeroN<0 | zeroN>0.5)) {
			stop("The elements of zeroN must be nonnegative and less than or equal 0.5.")
		}
	} else {
		zeroN<-1/(cbind(rep(R,S),Rp)*nstm)
		zeroN[zeroN==Inf]<-0
	}
	if(!is.null(start)) {
		if(length(start)!=(S*(R-1))) {
			stop("start must have the starting vector of the S(R-1) probabilities")
		}
		if(any( (start<0) | (start>1) )) {
			stop("start has improper starting values for some probabilities (outside the parameter space)")
		}
		if(any(kronecker(diag(S),rep(1,R-1))%*%start>rep(1,S))) {
			stop("The sum of the starting values for each subpopulation must not exceed 1")
		}
	} else {
		for(s in 1:S) {
			pp<-pst[[paste("st",s,".",1,sep="")]]
			if(any(pp==0)) {
				pp<-Nst[[paste("st",s,".",1,sep="")]]
				pp[pp==0]<-zero[s]
				start[((s-1)*(R-1)+1):(s*(R-1))]<-( pp/sum(pp) )[1:(R-1)]
			} else {
				start[((s-1)*(R-1)+1):(s*(R-1))]<-pp[1:(R-1)]
			}
		}
	}
	if(!is.numeric(maxit) || maxit<1) {
		stop("maxit must be >= 1")
	}
	if(!any(trace==c(0,1,2))) {
		stop("trace must be 0, 1 or 2")
	}
	if(!is.numeric(epsilon1) || epsilon1<=0) {
		stop("epsilon1 must be > 0")
	}
	if(!is.numeric(epsilon2) || epsilon2<=0 || epsilon2>=0.1) {
		stop("epsilon2 must be > 0 and < 0.1")
	}

	#Fisher information matrix i_1(\bar{theta})
	IF1<-function(thetab,missing,R,S,Tt,Zbst,Nst,nstm,nsmm) {
		if1<-matrix(0,S*(R-1),S*(R-1))
		for(s in 1:S) {
			for(tt in 1:Tt[s]) {
				if(tt==1) {
					zbst<-diag(R-1)
					rr<-R-1
				} else {
					zbst<-Zbst[[paste("st",s,".",tt,sep="")]]
					rr<-Rp[s,tt-1]-1
				}
				uns<-rep(1,rr)
				thetabst<-c(t(zbst) %*% thetab[((s-1)*(R-1)+1):(s*(R-1))])
				#\{\hat{alpha}_{st}\}
				if(missing=="MAR") {
					alphast<-c(solve(diag(c(thetabst,1-sum(thetabst))),tol=tole)%*%Nst[[paste("st",s,".",tt,sep="")]])/nsmm[s]
					if1[((s-1)*(R-1)+1):(s*(R-1)),((s-1)*(R-1)+1):(s*(R-1))]<-if1[((s-1)*(R-1)+1):(s*(R-1)),((s-1)*(R-1)+1):(s*(R-1))]+
					zbst%*%( diag(alphast[1:rr],nrow=rr)%*%solve(diag(thetabst,nrow=length(thetabst)),tol=tole) + (alphast[rr+1]/(1-sum(thetabst)))*uns%*%t(uns) )%*%t(zbst)
				}
				if(missing=="MCAR") {
					alphast<-nstm[s,tt]/nsmm[s]
					if1[((s-1)*(R-1)+1):(s*(R-1)),((s-1)*(R-1)+1):(s*(R-1))]<-if1[((s-1)*(R-1)+1):(s*(R-1)),((s-1)*(R-1)+1):(s*(R-1))]+
					alphast*zbst%*%( solve(diag(thetabst,nrow=length(thetabst)),tol=tole) + (1/(1-sum(thetabst)))*uns%*%t(uns) )%*%t(zbst)
				}
			}
			if1[((s-1)*(R-1)+1):(s*(R-1)),((s-1)*(R-1)+1):(s*(R-1))]<-if1[((s-1)*(R-1)+1):(s*(R-1)),((s-1)*(R-1)+1):(s*(R-1))]*nsmm[s]
		}
		return(if1)
	}

	#Likelihood ratio statistic of MCAR given MAR
	QvMCAR.MAR<-function(thetab,R,S,Zs,Ns,ps,bs,Bs) {
		Qv<-0
		for(s in 1:S) {
			Qv<-Qv-2*c( t(Ns[[paste("s",s,sep="")]][Ns[[paste("s",s,sep="")]]!=0]) %*%
				   (log(c(t(Zs[[paste("s",s,sep="")]])%*%(bs+Bs%*%thetab[((s-1)*(R-1)+1):(s*(R-1))]))[Ns[[paste("s",s,sep="")]]!=0])
				    -log(ps[[paste("s",s,sep="")]][Ns[[paste("s",s,sep="")]]!=0]) ) )
		}
		return(Qv)
	}

	conv<-FALSE
	it<-0
	thetab<-start
	theta<-b+c(B%*%thetab)
	QvMCAR<-QvMCAR.MAR(thetab,R,S,Zs,Ns,ps,bs,Bs)
	if(trace==1) {
		cat("Iteration - Likelihood ratio statistic of MCAR given MAR\n")
		cat(it," - ",QvMCAR,"\n")
	}
	if(trace==2) {
		cat("Iteration ",it," - LRS(MCAR|MAR) ",QvMCAR," - Estimates:\n",sep="")
		print.default(format(round(matrix(theta,S,R,byrow=T),digits=digits),digits=digits),quote=FALSE,print.gap=2)
		cat("\n")
	}
	if(any(method==c("FS-MCAR","NR/FS-MAR"))) {
		#Score vector S_1(\bar{theta})
		S1<-function(thetab,R,S,Tt,Zbst,pbst,nstm) {
			s1<-numeric(S*(R-1))
			for(s in 1:S) {
				for(tt in 1:Tt[s]) {
					if(tt==1) {
						zbst<-diag(R-1)
					} else {
						zbst<-Zbst[[paste("st",s,".",tt,sep="")]]
					}
					thetabst<-c(t(zbst) %*% thetab[((s-1)*(R-1)+1):(s*(R-1))])
					vthetabst<-(diag(thetabst,nrow=length(thetabst))-thetabst%*%t(thetabst))/nstm[s,tt]
					s1[((s-1)*(R-1)+1):(s*(R-1))]<-s1[((s-1)*(R-1)+1):(s*(R-1))]+
					c( zbst %*% solve(vthetabst,tol=tole) %*% (pbst[[paste("st",s,".",tt,sep="")]]-thetabst) )
				}
			}
			return(s1)
		}
	}
	while((conv==FALSE) & (it<maxit)) {
		it<-it+1
		thetaba<-thetab
		QvMCARa<-QvMCAR
		if(any(method==c("FS-MCAR","NR/FS-MAR"))) {
			if(method=="NR/FS-MAR") {
				IF1it<-IF1(thetab,missing="MAR",R,S,Tt,Zbst,Nst,nstm,nsmm)
			} else {
				IF1it<-IF1(thetab,missing="MCAR",R,S,Tt,Zbst,Nst,nstm,nsmm)
			}
			thetab<-thetaba+c( solve(IF1it,tol=tole) %*% S1(thetab,R,S,Tt,Zbst,pbst,nstm) )
			theta<-b+c(B%*%thetab)
		}
		if(method=="EM") {
			pp<-numeric(S*R)
			for(s in 1:S) {
				pp[((s-1)*R+1):(s*R)]<-Nst[[paste("st",s,".",1,sep="")]]
				if(Tt[s]>1) {
					for(tt in 2:Tt[s]) {
						zst<-Zst[[paste("st",s,".",tt,sep="")]]
						pp[((s-1)*R+1):(s*R)]<-pp[((s-1)*R+1):(s*R)]+
						diag(theta[((s-1)*R+1):(s*R)],nrow=R)%*%zst%*%
						solve(diag(c(t(zst)%*%theta[((s-1)*R+1):(s*R)]),nrow=Rp[s,tt-1]),tol=tole)%*%Nst[[paste("st",s,".",tt,sep="")]]
					}
				}
				pp[((s-1)*R+1):(s*R)]<-pp[((s-1)*R+1):(s*R)]/nsmm[s]
			}
			theta<-pp
			thetab<-c( kronecker(diag(S),cbind(diag(R-1),rep(0,R-1)))%*%theta )
		}
		if(any((theta<=0) | (theta>=1))) {
			if(any((theta<0) | (theta>1))) {
				stop(paste("Any of the estimated probabilities obtained by the iterative process are outside the parameter space.\nTry another iterative process or starting values. (iteration ",it,")",sep=""))
			} else {
				if(min(theta)==0) {
					warning(paste("Any of the estimated probabilities obtained by the iterative process are in the\nborderline of the parameter space (are equal 0). (iteration ",it,")",sep=""))
				}
				if(max(theta)==1) {
					warning(paste("Any of the estimated probabilities obtained by the iterative process are in the\nborderline of the parameter space (are equal 1). (iteration ",it,")",sep=""))
				}
			}
		}
		QvMCAR<-QvMCAR.MAR(thetab,R,S,Zs,Ns,ps,bs,Bs)
		if(trace==1) {
			cat(it," - ",QvMCAR,"\n")
		}
		if(trace==2) {
			cat("Iteration ",it," - LRS(MCAR|MAR) ",QvMCAR," - Estimates:\n",sep="")
			print.default(format(round(matrix(theta,S,R,byrow=T),digits=digits),digits=digits),quote=FALSE,print.gap=2)
			cat("\n")
		}
		if(abs(QvMCARa-QvMCAR)<epsilon1 && max(abs(thetab-thetaba))<epsilon2) {
			conv<-TRUE
		}
	}
	if(conv==FALSE) {
		warning(paste("\nThe iterative process has not attained the convergence criterions with ",maxit," iterations.\n",sep=""))
	}
	Vthetab<-solve(IF1(thetab,missing,R,S,Tt,Zbst,Nst,nstm,nsmm),tol=tole)
	Vtheta<-B%*%Vthetab%*%t(B)
	QpMCAR<-0 #Pearson statistic of MCAR given MAR
	QnMCAR<-0 #Neyman statistic of MCAR given MAR
	for(s in 1:S) {
		pp<-pst[[paste("st",s,".",1,sep="")]]
		if(any(pp==0)) {
			pp<-Nst[[paste("st",s,".",1,sep="")]]
			pp[pp==0]<-zeroN[s,1]
			pp<-pp/sum(pp)
		}
		if(Tt[s]>1) {
			for(tt in 2:Tt[s]) {
				pp1<-pst[[paste("st",s,".",tt,sep="")]]
				if(any(pp1==0)) {
					pp1<-Nst[[paste("st",s,".",tt,sep="")]]
					pp1[pp1==0]<-zeroN[s,tt]
					pp1<-pp1/sum(pp1)
				}
				pp<-c(pp,pp1)
			}
		}
		QpMCAR<-QpMCAR+c( t(ps[[paste("s",s,sep="")]]-c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)]))%*%
				 (diag(Nsm[[paste("s",s,sep="")]])%*%solve(diag(c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)])),tol=tole))%*%
		(ps[[paste("s",s,sep="")]]-c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)])) )
		QnMCAR<-QnMCAR+c( t(ps[[paste("s",s,sep="")]]-c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)]))%*%
				 (diag(Nsm[[paste("s",s,sep="")]])%*%solve(diag(pp),tol=tole))%*%
		(ps[[paste("s",s,sep="")]]-c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)])) )
	}
	glMCAR<-S+sum(l)-sum(Tt)
	alphast<-Nst
	yst<-Nst
	for(s in 1:S) {
		for(tt in 1:Tt[s]) {
			if(missing=="MAR") {
				if(tt==1) {
					zst<-diag(R)
				} else {
					zst<-Zst[[paste("st",s,".",tt,sep="")]]
				}
				thetast<-c(t(zst) %*% theta[((s-1)*R+1):(s*R)])
				alphast[[paste("st",s,".",tt,sep="")]]<-c(solve(diag(thetast),tol=tole)%*%Nst[[paste("st",s,".",tt,sep="")]])/nsmm[s]
				yst[[paste("st",s,".",tt,sep="")]]<-c( nsmm[s]*diag(theta[((s-1)*R+1):(s*R)])%*%zst%*%alphast[[paste("st",s,".",tt,sep="")]] )
			}
			if(missing=="MCAR") {
				alphast[[paste("st",s,".",tt,sep="")]]<-nstm[s,tt]/nsmm[s]
				yst[[paste("st",s,".",tt,sep="")]]<-nstm[s,tt]*theta[((s-1)*R+1):(s*R)]
			}
		}
	}
	res<-list(catdataobj=catdataobj,missing=missing,method=method,start=start,R=R,S=S,miss=miss,Rp=Rp,l=l,Tt=Tt,
		  Zs=Zs,Zst=Zst,Zbs=Zbs,Zbst=Zbst,Ns=Ns,Nst=Nst,Nsm=Nsm,nstm=nstm,nsmm=nsmm,ps=ps,pst=pst,pbst=pbst,bs=bs,Bs=Bs,b=b,B=B,
		  zero=zero,maxit=maxit,epsilon1=epsilon1,epsilon2=epsilon2,zeroN=zeroN,call=match.call(),conv=conv,it=it,
		  theta=theta,Vtheta=Vtheta,thetab=thetab,Vthetab=Vthetab,QvMCAR=QvMCAR,QpMCAR=QpMCAR,QnMCAR=QnMCAR,glMCAR=glMCAR,
		  alphast=alphast,yst=yst)
	class(res)<-"satMarML"
	res
}



print.satMarML<-function(x, digits=max(3,getOption("digits")-3), ...) {
	cat("\nCall: ",deparse(x$call),"\n",sep="")

	cat("\nS=",x$S," subpopulations  x  R=",x$R," response categories\n",sep="")

	pp<-matrix(0,x$S,x$R)
	epp<-matrix(0,x$S,x$R)
	for(s in 1:x$S) {
		pp[s,]<-x$theta[((s-1)*x$R+1):(s*x$R)]
		epp[s,]<-sqrt(diag(x$Vtheta[((s-1)*x$R+1):(s*x$R),((s-1)*x$R+1):(s*x$R)]))
	}

	cat("\nMaximum likelihood estimates of the probabilities:\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	cat("\nStandard errors (",x$missing,"):\n",sep="")
	print.default(format(round(epp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	pp<-matrix(0,3,2)
	pp[1,]<-c(x$QvMCAR,1-pchisq(x$QvMCAR,x$glMCAR))
	pp[2,]<-c(x$QpMCAR,1-pchisq(x$QpMCAR,x$glMCAR))
	pp[3,]<-c(x$QnMCAR,1-pchisq(x$QnMCAR,x$glMCAR))
	rownames(pp)<-c("Likelihood ratio","Pearson","Neyman")
	colnames(pp)<-c("statistic","p-value")

	cat("\nGoodness of fit statistics of MCAR given MAR assumption (d.f.=",x$glMCAR,")\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)

	cat("\n")
}


summary.satMarML<-function(object, digits=max(3,getOption("digits")-3), ...) {

	cat("\nCall: ",deparse(object$call),"\n\n",sep="")
	cat("\nS=",object$S," subpopulations  x  R=",object$R," response categories\n",sep="")

	pp<-matrix(0,object$S,object$R)
	epp<-matrix(0,object$S,object$R)
	for(s in 1:object$S) {
		pp[s,]<-object$theta[((s-1)*object$R+1):(s*object$R)]
		epp[s,]<-sqrt(diag(object$Vtheta[((s-1)*object$R+1):(s*object$R),((s-1)*object$R+1):(s*object$R)]))
	}

	cat("\n\nMaximum likelihood estimates of the probabilities:\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	cat("\nStandard errors (",object$missing,"):\n",sep="")
	print.default(format(round(epp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	if(object$conv==FALSE) {
		cat("\n\n",object$method," has NOT attained the convergence criterion with ",object$maxit," iterations.\n\n",sep="")
	} else {
		cat("\n\n",object$method," attained the convergence criterion in ",object$it," iterations.\n\n",sep="")
	}

	pp<-matrix(0,3,2)
	pp[1,]<-c(object$QvMCAR,1-pchisq(object$QvMCAR,object$glMCAR))
	pp[2,]<-c(object$QpMCAR,1-pchisq(object$QpMCAR,object$glMCAR))
	pp[3,]<-c(object$QnMCAR,1-pchisq(object$QnMCAR,object$glMCAR))
	rownames(pp)<-c("Likelihood ratio","Pearson","Neyman")
	colnames(pp)<-c("statistic","p-value")

	cat("\nGoodness of fit statistics of MCAR given MAR assumption (d.f.=",object$glMCAR,"):\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)

	cat("\n\nAugmented estimated frequencies under ",object$missing,":\n",sep="")
	for(s in 1:object$S) {
		if(object$S>1) {
			cat("\nSubpopulation ",s,"\n",sep="")
		}
		pp<-matrix(0,object$Tt[s],object$R)
		pp[1,]<-object$yst[[paste("st",s,".",1,sep="")]]
		if(object$Tt[s]>1) {
			for(tt in 2:object$Tt[s]) {
				pp[tt,]<-object$yst[[paste("st",s,".",tt,sep="")]]
			}
		}
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)
	}

	cat("\n")
}





#Fit saturated models for the marginal probabilities of categorization under MCAR by weighted least squares
satMcarWLS<-function(catdataobj,zeroN=NULL) {
	if(class(catdataobj)!="readCatdata") {
		stop("x is not of the class readCatdata")
	}
	tole<-1e-323 #tol argument for solve function

	R<-catdataobj$R
	S<-catdataobj$S
	miss<-catdataobj$miss
	Rp<-catdataobj$Rp
	l<-catdataobj$l
	Tt<-catdataobj$Tt
	Zbs<-catdataobj$Zbs
	Nst<-catdataobj$Nst
	nstm<-catdataobj$nstm
	pst<-catdataobj$pst
	pbs<-catdataobj$pbs
	pbst<-catdataobj$pbst
	Vpbst<-catdataobj$Vpbst
	b<-catdataobj$b
	B<-catdataobj$B

	if(miss==FALSE) {
		stop("This function is just to categorical data with missings.\nThe proportions in your object are already ML estimates of the probabilities.")
	}
	if(!is.null(zeroN)) {
		if(!is.matrix(zeroN) || (nrow(zeroN)!=S) || (ncol(zeroN)!=max(Tt))) {
			if(is.vector(zeroN,mode="numeric")) {
				zeroN<-matrix(zeroN[1],S,max(Tt))
			} else {
				stop("zeroN must be a number or a matrix with dimensions S x (max(Tt))")
			}
		}
		if(any(zeroN<0 | zeroN>0.5)) {
			stop("The elements of zeroN must be nonnegative and less than or equal 0.5.")
		}
	} else {
		zeroN<-1/(cbind(rep(R,S),Rp)*nstm)
		zeroN[zeroN==Inf]<-0
	}
	thetab<-numeric(S*(R-1))
	Vthetab<-matrix(0,S*(R-1),S*(R-1))
	VpbsI<-Vthetabs<-pbs
	for(s in 1:S) {
		VpbsI[[paste("s",s,sep="")]]<-matrix(0,R+l[s]-Tt[s],R+l[s]-Tt[s])
		i1<-1
		for(tt in 1:Tt[s]) {
			if(tt==1) {
				rr<-R-1
			} else {
				rr<-Rp[s,tt-1]-1
			}
			pp<-pst[[paste("st",s,".",tt,sep="")]]
			if(any(pp==0)) {
				pp<-Nst[[paste("st",s,".",tt,sep="")]]
				pp[pp==0]<-zeroN[s,tt]
				pp<-pp/sum(pp)
				pp<-pp[1:rr]
				VpbsI[[paste("s",s,sep="")]][i1:(i1+rr-1),i1:(i1+rr-1)]<-(diag(pp,nrow=length(pp))-pp%*%t(pp))/nstm[s,tt]
			} else {
				VpbsI[[paste("s",s,sep="")]][i1:(i1+rr-1),i1:(i1+rr-1)]<-Vpbst[[paste("st",s,".",tt,sep="")]]
			}
			i1<-i1+rr
		}
		VpbsI[[paste("s",s,sep="")]]<-solve(VpbsI[[paste("s",s,sep="")]],tol=tole)
		Vthetab[((s-1)*(R-1)+1):(s*(R-1)),((s-1)*(R-1)+1):(s*(R-1))]<-
			solve(Zbs[[paste("s",s,sep="")]]%*%VpbsI[[paste("s",s,sep="")]]%*%t(Zbs[[paste("s",s,sep="")]]),tol=tole)
		thetab[((s-1)*(R-1)+1):(s*(R-1))]<-c(Vthetab[((s-1)*(R-1)+1):(s*(R-1)),((s-1)*(R-1)+1):(s*(R-1))]%*%
						     Zbs[[paste("s",s,sep="")]]%*%VpbsI[[paste("s",s,sep="")]]%*%pbs[[paste("s",s,sep="")]])
	}
	theta<-b+c(B%*%thetab)
	if(any((theta<0) | (theta>1))) {
		warning("Any of the estimated probabilities are outside the parameter space.\n")
	}
	Vtheta<-B%*%Vthetab%*%t(B)
	QnMCAR<-0 #Neyman statistic of MCAR
	for(s in 1:S) {
		QnMCAR<-QnMCAR+c( t(pbs[[paste("s",s,sep="")]]-c(t(Zbs[[paste("s",s,sep="")]])%*%thetab[((s-1)*(R-1)+1):(s*(R-1))]))%*%
				 VpbsI[[paste("s",s,sep="")]]%*%(pbs[[paste("s",s,sep="")]]-c(t(Zbs[[paste("s",s,sep="")]])%*%thetab[((s-1)*(R-1)+1):(s*(R-1))])) )
	}
	glMCAR<-S+sum(l)-sum(Tt)
	yst<-catdataobj$Nst
	for(s in 1:S) {
		for(tt in 1:Tt[s]) {
			yst[[paste("st",s,".",tt,sep="")]]<-nstm[s,tt]*theta[((s-1)*R+1):(s*R)]
		}
	}
	res<-list(catdataobj=catdataobj,R=R,S=S,miss=miss,Rp=Rp,l=l,Tt=Tt,Zbs=Zbs,nstm=nstm,pbs=pbs,pbst=pbst,Vpbst,b=b,B=B,
		  zeroN=zeroN,call=match.call(),theta=theta,Vtheta=Vtheta,thetab=thetab,Vthetab=Vthetab,
		  QnMCAR=QnMCAR,glMCAR=glMCAR,yst=yst)
	class(res)<-"satMcarWLS"
	res
}

print.satMcarWLS<-function(x, digits=max(3,getOption("digits")-3), ...) {
	cat("\nCall: ",deparse(x$call),"\n",sep="")

	cat("\nS=",x$S," subpopulations  x  R=",x$R," response categories\n",sep="")

	pp<-matrix(0,x$S,x$R)
	epp<-matrix(0,x$S,x$R)
	for(s in 1:x$S) {
		pp[s,]<-x$theta[((s-1)*x$R+1):(s*x$R)]
		epp[s,]<-sqrt(diag(x$Vtheta[((s-1)*x$R+1):(s*x$R),((s-1)*x$R+1):(s*x$R)]))
	}

	cat("\nWeighted least squares estimates of the probabilities:\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	cat("\nStandard errors:\n",sep="")
	print.default(format(round(epp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	cat("\nNeyman goodness of fit statistic of MCAR (d.f.=",x$glMCAR,"): ",
	    round(x$QnMCAR,digits=digits)," (p-value=",
	    format(round(1-pchisq(x$QnMCAR,x$glMCAR),digits=digits),digits=digits),")\n",sep="")

	cat("\n")
}


summary.satMcarWLS<-function(object, digits=max(3,getOption("digits")-3), ...) {

	cat("\nCall: ",deparse(object$call),"\n\n",sep="")
	cat("\nS=",object$S," subpopulations  x  R=",object$R," response categories\n",sep="")

	pp<-matrix(0,object$S,object$R)
	epp<-matrix(0,object$S,object$R)
	for(s in 1:object$S) {
		pp[s,]<-object$theta[((s-1)*object$R+1):(s*object$R)]
		epp[s,]<-sqrt(diag(object$Vtheta[((s-1)*object$R+1):(s*object$R),((s-1)*object$R+1):(s*object$R)]))
	}

	cat("\n\nWeighted least squares estimates of the probabilities:\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	cat("\nStandard errors:\n",sep="")
	print.default(format(round(epp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	cat("\n\nNeyman goodness of fit statistic of MCAR (d.f.=",object$glMCAR,"): ",
	    round(object$QnMCAR,digits=digits)," (p-value=",
	    format(round(1-pchisq(object$QnMCAR,object$glMCAR),digits=digits),digits=digits),")\n",sep="")

	cat("\n\nAugmented estimated frequencies under MCAR:\n",sep="")
	for(s in 1:object$S) {
		if(object$S>1) {
			cat("\nSubpopulation ",s,"\n",sep="")
		}
		pp<-matrix(0,object$Tt[s],object$R)
		pp[1,]<-object$yst[[paste("st",s,".",1,sep="")]]
		if(object$Tt[s]>1) {
			for(tt in 2:object$Tt[s]) {
				pp[tt,]<-object$yst[[paste("st",s,".",tt,sep="")]]
			}
		}
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)
	}
	cat("\n")
}





#Fit (general) linear models by maximum likelihood to the marginal probabilities of categorization for complete data or under MAR/MCAR
linML<-function(obj,A=NULL,X=NULL,U=NULL,start=NULL,maxit=100,trace=0,epsilon1=1e-6,epsilon2=1e-6,zeroN=NULL, digits=max(3,getOption("digits")-3)) {
	if(class(obj)!="readCatdata" && class(obj)!="satMarML") {
		stop("obj is not of the class readCatdata or satMarML")
	}
	if(class(obj)=="readCatdata" && obj$miss==TRUE) {
		stop("With missing data, obj must be of class satMarML")
	}
	tole<-1e-323 #tol argument for solve function

	R<-obj$R
	S<-obj$S
	miss<-obj$miss
	Rp<-obj$Rp
	l<-obj$l
	Tt<-obj$Tt
	Zs<-obj$Zs
	Zbs<-obj$Zbs
	Ns<-obj$Ns
	Nst<-obj$Nst
	Nsm<-obj$Nsm
	nstm<-obj$nstm
	nsmm<-obj$nsmm
	ps<-obj$ps
	pst<-obj$pst
	pbst<-obj$pbst
	bs<-obj$bs
	Bs<-obj$Bs
	b<-obj$b
	B<-obj$B
	theta<-obj$theta
	Vtheta<-obj$Vtheta
	thetab<-obj$thetab
	Vthetab<-obj$Vthetab

	if(miss==TRUE) {
		Zst<-obj$Zst
		Zbst<-obj$Zbst
		missing<-obj$missing
		QvMCAR<-obj$QvMCAR
		alphast<-obj$alphast
	} else {
		missing<-"MCAR"
		alphast<-Nst
		for(s in 1:S) {
			alphast[[s]]<-1
		}
	}

	restr<-kronecker(diag(S),rep(1,R))
	if(is.null(A)) {
		A<-kronecker(diag(S),cbind(diag(R-1),rep(0,R-1)))
		u<-S*(R-1)
	} else {
		if(!is.matrix(A)) {
			A<-t(A)
		}
		if(ncol(A)!=S*R) {
			stop(paste("The number of columns of A must be equal to S*R=",S*R,".",sep=""))
		}
		u<-nrow(A)
		if(u>S*(R-1)) {
			stop(paste("The number of rows of A must be less than or equal to S*(R-1)=",S*(R-1),".",sep=""))
		}
		if(qr(A)$rank!=u) {
			stop("Rank of A and its number of rows must be the same.")
		}
		if(qr(cbind(t(A),restr))$rank!=(u+S)) {
			stop("The lines of A are not linearly independent from the columns\nof the natural restrictions (kronecker(diag(S),rep(1,R))).")
		}
	}
	if(is.null(X) && is.null(U)) {
		stop("Specify X for freedom equation formulation or U for constraint formulation.")
	}
	if((!is.null(X)) && (!is.null(U))) {
		stop("Specify just X or U, not both.")
	}

	#Obtém a matriz U (ou X) a partir da X (ou U), via decomposicao valor singular de uma matriz Q
	#ortogonal a X (ou U). A matriz Q, ortogonal a X, e da forma: Q = Ic - X*inv(X'X)*X'
	#Q possui posto c-r, onde r e o posto de X.
	COMPLEMENT<-function(A) {
		K<-nrow(A)
		J<-ncol(A)
		if(J!=K) { #A nao e' quadrada de posto completo */
			if(J > K) {
				A<-t(A)
				K<-nrow(A)
				J<-ncol(A)
			}
			Q<-diag(K)-A%*%solve(t(A)%*%A,tol=tole)%*%t(A)
			Q.svd<-svd(Q) #decomposicao do valor singular, Q=Q.svd$u %*% diag(Q.svd$d) %*% t(Q.svd$v) ou Q=t(Q.svd$v %*% (t(Q.svd$u) * Q.svd$d))
			#base ortogonal
			B1<-Q.svd$d
			B2<-Q.svd$u
			B3<-Q.svd$v
			CAUX<-matrix(0, nrow=nrow(B2), ncol=1)
			for(i in 1:length(B1)) {
				if(round(B1[i],10)==1) {
					CAUX<-cbind(CAUX,B2[,i])
				}
			}
			if(J < K) {
				CAUX<-CAUX[,2:ncol(CAUX)]
			}
		} else {
			CAUX<-0
		}
		return(CAUX)
	}

	if(is.null(U)) {
		if(!is.matrix(X)) {
			X<-as.matrix(X)
		}
		if(nrow(X)!=u) {
			stop("X and A must have the same number of rows.")
		}
		p<-ncol(X)
		if(p>u) {
			stop("X can't have more columns than rows.")
		}
		if(qr(X)$rank!=p) {
			stop("Rank of X and its number of columns must be the same.")
		}
		U<-t(COMPLEMENT(X))
		form<-"freedom equation"
		if(!is.matrix(U)) {
			U<-t(U)
		}
	} else {
		if(!is.matrix(U)) {
			U<-t(U)
		}
		if(ncol(U)!=u) {
			stop("The number of columns of U have to be the same of the number of rows of A.")
		}
		p<-u-nrow(U)
		if(p==0) {
			stop("The number of rows and columns of U must not be the same.")
		}
		if(qr(U)$rank!=(u-p)) {
			stop("Rank of U and its number of rows must be the same.")
		}
		X<-COMPLEMENT(U)
		form<-"constraint"
		if(!is.matrix(X)) {
			X<-as.matrix(X)
		}
	}

	if(u<S*(R-1)) {
		A0<-t(COMPLEMENT(rbind(A,t(restr))))
		Wpart<-kronecker(diag(S),cbind(diag(R-1),rep(0,R-1)))%*%solve(rbind(A,t(restr),A0),tol=tole)
		W<-Wpart%*%rbind(cbind(X,matrix(0,u,S*(R-1)-u)),matrix(0,S,S*(R-1)-u+p),cbind(matrix(0,S*(R-1)-u,p),diag(S*(R-1)-u)))
	} else {
		Wpart<-kronecker(diag(S),cbind(diag(R-1),rep(0,R-1)))%*%solve(rbind(A,t(restr)),tol=tole)
		W<-Wpart%*%rbind(X,matrix(0,S,p))
	}

	if(!is.null(zeroN)) {
		if(!is.matrix(zeroN) || (nrow(zeroN)!=S) || (ncol(zeroN)!=max(Tt))) {
			if(is.vector(zeroN,mode="numeric")) {
				zeroN<-matrix(zeroN[1],S,max(Tt))
			} else {
				stop("zeroN must be a number or a matrix with dimensions S x (max(Tt))")
			}
		}
		if(any(zeroN<0 | zeroN>0.5)) {
			stop("The elements of zeroN must be nonnegative and less than or equal 0.5.")
		}
	} else {
		if(miss==TRUE) {
			zeroN<-1/(cbind(rep(R,S),Rp)*nstm)
			zeroN[zeroN==Inf]<-0
		} else {
			zeroN<-1/(rep(R,S)*nstm)
		}
	}
	if(!is.numeric(maxit) || maxit<1) {
		stop("maxit must be >= 1")
	}
	if(!any(trace==c(0,1,2))) {
		stop("trace must be 0, 1 or 2")
	}
	if(!is.numeric(epsilon1) || epsilon1<=0) {
		stop("epsilon1 must be > 0")
	}
	if(!is.numeric(epsilon2) || epsilon2<=0 || epsilon2>=0.1) {
		stop("epsilon2 must be > 0 and < 0.1")
	}

	#Score vector S_1(\bar{theta})
	S1<-function(thetab,R,S,Tt,Zbst,pbst,nstm) {
		s1<-numeric(S*(R-1))
		for(s in 1:S) {
			for(tt in 1:Tt[s]) {
				if(tt==1) {
					zbst<-diag(R-1)
				} else {
					zbst<-Zbst[[paste("st",s,".",tt,sep="")]]
				}
				thetabst<-c(t(zbst) %*% thetab[((s-1)*(R-1)+1):(s*(R-1))])
				vthetabst<-(diag(thetabst,nrow=length(thetabst))-thetabst%*%t(thetabst))/nstm[s,tt]
				s1[((s-1)*(R-1)+1):(s*(R-1))]<-s1[((s-1)*(R-1)+1):(s*(R-1))]+
				c( zbst %*% solve(vthetabst,tol=tole) %*% (pbst[[paste("st",s,".",tt,sep="")]]-thetabst) )
			}
		}
		return(s1)
	}

	#Fisher information matrix i_1(\bar{theta})
	#>>>>>>>>>> Esta informação de Fisher é diferente da existente em satMarML !!!
	#Aquela estima iterativamente alphast e esta utiliza sempre o alphast obtido de theta sem restrições
	IF1<-function(thetab,missing,R,S,Tt,Zbst,Nst,nstm,nsmm,alphast) {
		if1<-matrix(0,S*(R-1),S*(R-1))
		for(s in 1:S) {
			for(tt in 1:Tt[s]) {
				if(tt==1) {
					zbst<-diag(R-1)
					rr<-R-1
				} else {
					zbst<-Zbst[[paste("st",s,".",tt,sep="")]]
					rr<-Rp[s,tt-1]-1
				}
				uns<-rep(1,rr)
				thetabst<-c(t(zbst) %*% thetab[((s-1)*(R-1)+1):(s*(R-1))])
				#\{\hat{alpha}_{st}\}
				if(missing=="MAR") {
					alpha<-alphast[[paste("st",s,".",tt,sep="")]]
					if1[((s-1)*(R-1)+1):(s*(R-1)),((s-1)*(R-1)+1):(s*(R-1))]<-if1[((s-1)*(R-1)+1):(s*(R-1)),((s-1)*(R-1)+1):(s*(R-1))]+
					zbst%*%( diag(alpha[1:rr],nrow=rr)%*%solve(diag(thetabst,nrow=length(thetabst)),tol=tole) + (alpha[rr+1]/(1-sum(thetabst)))*uns%*%t(uns) )%*%t(zbst)
				}
				if(missing=="MCAR") {
					alpha<-alphast[[paste("st",s,".",tt,sep="")]]
					if1[((s-1)*(R-1)+1):(s*(R-1)),((s-1)*(R-1)+1):(s*(R-1))]<-if1[((s-1)*(R-1)+1):(s*(R-1)),((s-1)*(R-1)+1):(s*(R-1))]+
					alpha*zbst%*%( solve(diag(thetabst,nrow=length(thetabst)),tol=tole) + (1/(1-sum(thetabst)))*uns%*%t(uns) )%*%t(zbst)
				}
			}
			if1[((s-1)*(R-1)+1):(s*(R-1)),((s-1)*(R-1)+1):(s*(R-1))]<-if1[((s-1)*(R-1)+1):(s*(R-1)),((s-1)*(R-1)+1):(s*(R-1))]*nsmm[s]
		}
		return(if1)
	}

	#Likelihood ratio statistic of H given MAR or any missing mechanism contained in MAR
	Qv<-function(thetaH,theta,R,S,Zs,Ns) {
		Qv<-0
		for(s in 1:S) {
			Qv<-Qv-2*c( t(Ns[[paste("s",s,sep="")]][Ns[[paste("s",s,sep="")]]!=0]) %*%
				   (log(c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)])[Ns[[paste("s",s,sep="")]]!=0])
				    -log(c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)])[Ns[[paste("s",s,sep="")]]!=0]) ) )
		}
		return(Qv)
	}
	conv<-FALSE
	it<-0
	if(is.null(start)) {
		theta0<-theta
		Vtheta0<-Vtheta
		if(miss==FALSE) {
			for(s in 1:S) {
				if(any(theta[((s-1)*R+1):(s*R)]==0)) {
					theta0[((s-1)*R+1):(s*R)]<-Nst[[paste("st",s,".",1,sep="")]]
					theta0[((s-1)*R+1):(s*R)][theta0[((s-1)*R+1):(s*R)]==0]<-zeroN[s,1]
					theta0[((s-1)*R+1):(s*R)]<-theta0[((s-1)*R+1):(s*R)]/sum(theta0[((s-1)*R+1):(s*R)])
					Vtheta0[((s-1)*R+1):(s*R),((s-1)*R+1):(s*R)]<-(diag(theta0[((s-1)*R+1):(s*R)])-theta0[((s-1)*R+1):(s*R)]%*%t(theta0[((s-1)*R+1):(s*R)]))/nsmm[s]
				}
			}
		}
		if(u<S*(R-1)) {
			AA<-rbind(A,A0)
			XX<-rbind(cbind(X,matrix(0,u,S*(R-1)-u)),cbind(matrix(0,S*(R-1)-u,p),diag(S*(R-1)-u)))
		} else {
			AA<-A
			XX<-X
		}
		beta<-c( solve(t(XX)%*%solve(AA%*%Vtheta0%*%t(AA),tol=tole)%*%XX,tol=tole) %*% t(XX)%*%solve(AA%*%Vtheta0%*%t(AA),tol=tole)%*%AA%*%theta0 )
	} else {
		if(!is.vector(start,mode="numeric")) {
			stop("start must be a numeric vector")
		}
		if(u==S*(R-1)) {
			if(form=="constraint") {
				stop("start may be used just with the freedom equation formulation")
			}
			if(length(start)!=p) {
				stop("start must have the same length of the number of parameters")
			}
			beta<-start
		} else {
			stop("start may be supplied just if the number of rows of A is S(R-1)")
		}
	}
	if(u<S*(R-1)){ thetabH<-Wpart%*%c(X%*%beta[1:p],rep(1,S),beta[(p+1):(S*(R-1)-u+p)]) } else { thetabH<-Wpart%*%c(X%*%beta,rep(1,S)) }
	thetaH<-b+c(B%*%thetabH)
	QvH<-Qv(thetaH,theta,R,S,Zs,Ns)
	if(trace==1) {
		cat("Iteration - Likelihood ratio statistic of the LM\n")
		cat(it," - ",QvH,"\n")
	}
	if(trace==2) {
		cat("Iteration ",it," - LRS(LM) ",QvH," - Parameter estimates:\n",sep="")
		print.default(format(round(beta[1:p],digits=digits),digits=digits),quote=FALSE,print.gap=2)
		cat("\n")
	}
	while((conv==FALSE) & (it<maxit)) {
		it<-it+1
		betaa<-beta
		QvHa<-QvH
		beta<-betaa+c( solve(t(W)%*%IF1(thetabH,missing,R,S,Tt,Zbst,Nst,nstm,nsmm,alphast)%*%W,tol=tole) %*%
			      t(W) %*% S1(thetabH,R,S,Tt,Zbst,pbst,nstm) )
		if(u<S*(R-1)){ thetabH<-Wpart%*%c(X%*%beta[1:p],rep(1,S),beta[(p+1):(S*(R-1)-u+p)]) } else { thetabH<-Wpart%*%c(X%*%beta,rep(1,S)) }
		thetaH<-b+c(B%*%thetabH)
		if(any((thetaH<=0) | (thetaH>=1))) {
			if(any((thetaH<0) | (thetaH>1))) {
				stop(paste("Any of the estimated probabilities obtained by the iterative\nprocess are outside the parameter space. (iteration ",it,")",sep=""))
			} else {
				if(min(thetaH)==0) {
					warning(paste("Any of the estimated probabilities obtained by the iterative process are in the\nborderline of the parameter space (are equal 0). (iteration ",it,")",sep=""))
				}
				if(max(thetaH)==1) {
					warning(paste("Any of the estimated probabilities obtained by the iterative process are in the\nborderline of the parameter space (are equal 1). (iteration ",it,")",sep=""))
				}
			}
		}
		QvH<-Qv(thetaH,theta,R,S,Zs,Ns)
		if(trace==1) {
			cat(it," - ",QvH,"\n")
		}
		if(trace==2) {
			cat("Iteration ",it," - LRS(LM) ",QvH," - Parameter estimates:\n",sep="")
			print.default(format(round(beta[1:p],digits=digits),digits=digits),quote=FALSE,print.gap=2)
			cat("\n")
		}
		if(abs(QvHa-QvH)<epsilon1 && max(abs(beta-betaa))<epsilon2) {
			conv<-TRUE
		}
	}
	if(conv==FALSE) {
		warning(paste("\nThe iterative process has not attained the convergence criterions with ",maxit," iterations.\n",sep=""))
	}
	Vbeta<-solve(t(W)%*%IF1(thetabH,missing,R,S,Tt,Zbst,Nst,nstm,nsmm,alphast)%*%W,tol=tole)
	VthetabH<-W%*%Vbeta%*%t(W)
	VthetaH<-B%*%VthetabH%*%t(B)
	beta<-beta[1:p]
	Vbeta<-Vbeta[1:p,1:p]
	Fu<-c(A%*%theta) #F=A*theta unrestricted
	FH<-c(X%*%beta)  #F=A*theta under LM
	VFu<-A%*%Vtheta%*%t(A)
	VFH<-X%*%Vbeta%*%t(X)
	if(miss==TRUE) {
		QvHMCAR<-QvH+QvMCAR #Likelihood ratio statistic of H and MCAR given MAR
		QpHMCAR<-0 #Pearson statistic of H and MCAR given MAR
		QnHMCAR<-0 #Neyman statistic of H and MCAR given MAR
	}
	QpH<-0 #Pearson statistic of H given missing mechanism (MAR or MCAR)
	QnH<-0 #Neyman statistic of H given missing mechanism (MAR or MCAR)
	for(s in 1:S) {
		if(miss==TRUE) {
			pp<-pst[[paste("st",s,".",1,sep="")]]
			if(any(pp==0)) {
				pp<-Nst[[paste("st",s,".",1,sep="")]]
				pp[pp==0]<-zeroN[s,1]
				pp<-pp/sum(pp)
			}
			if(missing!="MAR") {
				th<-c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)])
			}
		} else {
			th<-pst[[paste("st",s,".",1,sep="")]]
			if(any(th==0)) {
				th<-Nst[[paste("st",s,".",1,sep="")]]
				th[th==0]<-zeroN[s,1]
				th<-th/sum(th)
			}
		}
		if(Tt[s]>1) {
			for(tt in 2:Tt[s]) {
				pp1<-pst[[paste("st",s,".",tt,sep="")]]
				if(any(pp1==0)) {
					pp1<-Nst[[paste("st",s,".",tt,sep="")]]
					pp1[pp1==0]<-zeroN[s,tt]
					pp1<-pp1/sum(pp1)
				}
				pp<-c(pp,pp1)
			}
		}
		if(missing=="MAR") {
			QpH<-QpH+c( t(c(t(Zs[[paste("s",s,sep="")]])%*%(theta[((s-1)*R+1):(s*R)]-thetaH[((s-1)*R+1):(s*R)])))%*%
				   (diag(Ns[[paste("s",s,sep="")]])%*%solve(diag(c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)])),tol=tole)%*%solve(diag(c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)])),tol=tole))%*%
			(c(t(Zs[[paste("s",s,sep="")]])%*%(theta[((s-1)*R+1):(s*R)]-thetaH[((s-1)*R+1):(s*R)]))) )
			QnH<-QnH+c( t(rep(1,R+l[s])-c(solve(diag(c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)])),tol=tole)%*%t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)]))%*%
				   diag(Ns[[paste("s",s,sep="")]])%*%
			(rep(1,R+l[s])-c(solve(diag(c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)])),tol=tole)%*%t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)])) )
		} else {
			QpH<-QpH+c( t(c(t(Zs[[paste("s",s,sep="")]])%*%(theta[((s-1)*R+1):(s*R)]-thetaH[((s-1)*R+1):(s*R)])))%*%
				   (diag(Nsm[[paste("s",s,sep="")]])%*%solve(diag(c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)])),tol=tole))%*%
			(c(t(Zs[[paste("s",s,sep="")]])%*%(theta[((s-1)*R+1):(s*R)]-thetaH[((s-1)*R+1):(s*R)]))) )
			QnH<-QnH+c( t(c(t(Zs[[paste("s",s,sep="")]])%*%(theta[((s-1)*R+1):(s*R)]-thetaH[((s-1)*R+1):(s*R)])))%*%
				   (diag(Nsm[[paste("s",s,sep="")]])%*%solve(diag(th),tol=tole))%*%
			(c(t(Zs[[paste("s",s,sep="")]])%*%(theta[((s-1)*R+1):(s*R)]-thetaH[((s-1)*R+1):(s*R)]))) )
		}
		if(miss==TRUE) {
			QpHMCAR<-QpHMCAR+c( t(ps[[paste("s",s,sep="")]]-c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)]))%*%
					   (diag(Nsm[[paste("s",s,sep="")]])%*%solve(diag(c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)])),tol=tole))%*%
			(ps[[paste("s",s,sep="")]]-c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)])) )
			QnHMCAR<-QnHMCAR+c( t(ps[[paste("s",s,sep="")]]-c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)]))%*%
					   (diag(Nsm[[paste("s",s,sep="")]])%*%solve(diag(pp),tol=tole))%*%
			(ps[[paste("s",s,sep="")]]-c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)])) )
		}
	}
	glH<-u-p
	glHMCAR<-glH+S+sum(l)-sum(Tt)
	QwH<-0
	if(glH>0) {
		QwH<-c(t(U%*%Fu)%*%solve(U%*%VFu%*%t(U),tol=tole)%*%U%*%Fu)
	}
	ystH<-Nst
	for(s in 1:S) {
		for(tt in 1:Tt[s]) {
			if(missing=="MAR") {
				if(tt==1) {
					zst<-diag(R)
				} else {
					zst<-Zst[[paste("st",s,".",tt,sep="")]]
				}
				ystH[[paste("st",s,".",tt,sep="")]]<-c( nsmm[s]*diag(thetaH[((s-1)*R+1):(s*R)])%*%zst%*%alphast[[paste("st",s,".",tt,sep="")]] )
			}
			if(missing=="MCAR") {
				ystH[[paste("st",s,".",tt,sep="")]]<-nstm[s,tt]*thetaH[((s-1)*R+1):(s*R)]
			}
		}
	}
	if(miss==TRUE) {
		res<-list(obj=obj,missing=missing,R=R,S=S,miss=miss,Rp=Rp,l=l,Tt=Tt,
			  Zs=Zs,Zst=Zst,Zbs=Zbs,Zbst=Zbst,Ns=Ns,Nst=Nst,Nsm=Nsm,nstm=nstm,nsmm=nsmm,ps=ps,pst=pst,pbst=pbst,bs=bs,Bs=Bs,b=b,B=B,
			  theta=theta,Vtheta=Vtheta,thetab=thetab,Vthetab=Vthetab,alphast=alphast,maxit=maxit,epsilon1=epsilon1,
			  epsilon2=epsilon2,zeroN=zeroN,call=match.call(),A=A,X=X,U=U,form=form,p=p,u=u,conv=conv,it=it,
			  thetaH=thetaH,VthetaH=VthetaH,thetabH=thetabH,VthetabH=VthetabH,beta=beta,Vbeta=Vbeta,Fu=Fu,VFu=VFu,FH=FH,VFH=VFH,
			  QvH=QvH,QpH=QpH,QnH=QnH,QwH=QwH,glH=glH,QvHMCAR=QvHMCAR,QpHMCAR=QpHMCAR,QnHMCAR=QnHMCAR,glHMCAR=glHMCAR,ystH=ystH)
	} else {
		res<-list(obj=obj,R=R,S=S,miss=miss,Rp=Rp,l=l,Tt=Tt,
			  Zs=Zs,Zbs=Zbs,Ns=Ns,Nst=Nst,Nsm=Nsm,nstm=nstm,nsmm=nsmm,ps=ps,pst=pst,pbst=pbst,bs=bs,Bs=Bs,b=b,B=B,
			  theta=theta,Vtheta=Vtheta,thetab=thetab,Vthetab=Vthetab,maxit=maxit,epsilon1=epsilon1,
			  epsilon2=epsilon2,zeroN=zeroN,call=match.call(),A=A,X=X,U=U,form=form,p=p,u=u,conv=conv,it=it,
			  thetaH=thetaH,VthetaH=VthetaH,thetabH=thetabH,VthetabH=VthetabH,beta=beta,Vbeta=Vbeta,Fu=Fu,VFu=VFu,FH=FH,VFH=VFH,
			  QvH=QvH,QpH=QpH,QnH=QnH,QwH=QwH,glH=glH,ystH=ystH)
	}
	class(res)<-"linML"
	res
}



print.linML<-function(x, digits=max(3,getOption("digits")-3), ...) {
	cat("\nCall: ",deparse(x$call),"\n",sep="")

	if(x$form=="freedom equation") {
		pp<-matrix(0,x$p,4)
		pp[,1]<-x$beta
		pp[,2]<-sqrt(diag(as.matrix(x$Vbeta)))
		pp[,3]<-pp[,1]/pp[,2]
		pp[,4]<-1-pchisq((pp[,3])^2,1)
		colnames(pp)<-c("estimate","std.error","z-value","p-value")
		cat("\nMaximum likelihood estimates of the parameters of the linear model",sep="")
		if(x$miss==TRUE) {
			cat(" under ",x$missing,sep="")
		}
		cat(":\n",sep="")
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)
	}

	pp<-matrix(0,4,2)
	pp[1,]<-c(x$QvH,ifelse(x$glH>0,1-pchisq(x$QvH,x$glH),1))
	pp[2,]<-c(x$QpH,ifelse(x$glH>0,1-pchisq(x$QpH,x$glH),1))
	pp[3,]<-c(x$QnH,ifelse(x$glH>0,1-pchisq(x$QnH,x$glH),1))
	pp[4,]<-c(x$QwH,ifelse(x$glH>0,1-pchisq(x$QwH,x$glH),1))
	rownames(pp)<-c("Likelihood ratio","Pearson","Neyman","Wald")
	colnames(pp)<-c("statistic","p-value")

	cat("\nGoodness of fit of the linear model",sep="")
	if(x$miss==TRUE) {
		cat(" given ",x$missing,sep="")
	}
	cat(" (d.f.=",x$glH,"):\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)

	if(x$miss==TRUE) {
		pp<-matrix(0,3,2)
		pp[1,]<-c(x$QvHMCAR,1-pchisq(x$QvHMCAR,x$glHMCAR))
		pp[2,]<-c(x$QpHMCAR,1-pchisq(x$QpHMCAR,x$glHMCAR))
		pp[3,]<-c(x$QnHMCAR,1-pchisq(x$QnHMCAR,x$glHMCAR))
		rownames(pp)<-c("Likelihood ratio","Pearson","Neyman")
		colnames(pp)<-c("statistic","p-value")
		cat("\nGoodness of fit of the linear model and MCAR given MAR (d.f.=",x$glHMCAR,"):\n",sep="")
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)
	}

	cat("\n")
}


summary.linML<-function(object, digits=max(3,getOption("digits")-3), ...) {

	cat("\nCall: ",deparse(object$call),"\n\n",sep="")

	pp<-matrix(0,object$S,object$R)
	epp<-matrix(0,object$S,object$R)
	for(s in 1:object$S) {
		pp[s,]<-object$thetaH[((s-1)*object$R+1):(s*object$R)]
		epp[s,]<-sqrt(diag(object$VthetaH[((s-1)*object$R+1):(s*object$R),((s-1)*object$R+1):(s*object$R)]))
	}
	cat("\nMaximum likelihood estimates of the probabilities under the linear model (LM):\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	cat("\nStandard errors",sep="")
	if(object$miss==TRUE) {
		cat(" (",object$missing,")",sep="")
	}
	cat(":\n",sep="")
	print.default(format(round(epp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	pp<-matrix(0,object$u,4)
	pp[,1]<-object$Fu
	pp[,2]<-sqrt(diag(as.matrix(object$VFu)))
	pp[,3]<-object$FH
	pp[,4]<-sqrt(diag(as.matrix(object$VFH)))
	colnames(pp)<-c("observed","std.error","under the LM","std.error")
	cat("\n\nMaximum likelihood estimates of the linear functions specified by the matrix A:\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	if(object$form=="freedom equation") {
		pp<-matrix(0,object$p,4)
		pp[,1]<-object$beta
		pp[,2]<-sqrt(diag(as.matrix(object$Vbeta)))
		pp[,3]<-pp[,1]/pp[,2]
		pp[,4]<-1-pchisq((pp[,3])^2,1)
		colnames(pp)<-c("estimate","std.error","z-value","p-value")
		cat("\n\nMaximum likelihood estimates of the parameters of the linear model",sep="")
		if(object$miss==TRUE) {
			cat(" under ",object$missing,sep="")
		}
		cat(":\n",sep="")
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)
	}

	if(object$conv==FALSE) {
		cat("\n\nFisher scoring has NOT attained the convergence criterion with ",object$maxit," iterations.\n",sep="")
	} else {
		cat("\n\nFisher scoring attained the convergence criterion in ",object$it," iterations.\n",sep="")
	}

	pp<-matrix(0,4,2)
	pp[1,]<-c(object$QvH,ifelse(object$glH>0,1-pchisq(object$QvH,object$glH),1))
	pp[2,]<-c(object$QpH,ifelse(object$glH>0,1-pchisq(object$QpH,object$glH),1))
	pp[3,]<-c(object$QnH,ifelse(object$glH>0,1-pchisq(object$QnH,object$glH),1))
	pp[4,]<-c(object$QwH,ifelse(object$glH>0,1-pchisq(object$QwH,object$glH),1))
	rownames(pp)<-c("Likelihood ratio","Pearson","Neyman","Wald")
	colnames(pp)<-c("statistic","p-value")

	cat("\n\nGoodness of fit of the linear model",sep="")
	if(object$miss==TRUE) {
		cat(" given ",object$missing,sep="")
	}
	cat(" (d.f.=",object$glH,"):\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)

	if(object$miss==TRUE) {
		pp<-matrix(0,3,2)
		pp[1,]<-c(object$QvHMCAR,1-pchisq(object$QvHMCAR,object$glHMCAR))
		pp[2,]<-c(object$QpHMCAR,1-pchisq(object$QpHMCAR,object$glHMCAR))
		pp[3,]<-c(object$QnHMCAR,1-pchisq(object$QnHMCAR,object$glHMCAR))
		rownames(pp)<-c("Likelihood ratio","Pearson","Neyman")
		colnames(pp)<-c("statistic","p-value")
		cat("\n\nGoodness of fit of the linear model and MCAR given MAR (d.f.=",object$glHMCAR,"):\n",sep="")
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)
	}

	if(object$miss==TRUE) {
		cat("\n\nAugmented estimated frequencies under the linear model and ",object$missing,":\n",sep="")
		for(s in 1:object$S) {
			if(object$S>1) {
				cat("\nSubpopulation ",s,"\n",sep="")
			}
			pp<-matrix(0,object$Tt[s],object$R)
			pp[1,]<-object$yst[[paste("st",s,".",1,sep="")]]
			if(object$Tt[s]>1) {
				for(tt in 2:object$Tt[s]) {
					pp[tt,]<-object$yst[[paste("st",s,".",tt,sep="")]]
				}
			}
			print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)
		}
	} else {
		cat("\n\nEstimated frequencies under the linear model:\n",sep="")
		pp<-matrix(0,object$S,object$R)
		for(s in 1:object$S) {
			pp[s,]<-object$yst[[paste("st",s,".",1,sep="")]]
		}
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)
	}
	cat("\n")
}





#Fit log-linear models by maximum likelihood to the marginal probabilities of categorization for complete data or under MAR/MCAR
loglinML<-function(obj,A=NULL,X=NULL,U=NULL,XL=NULL,UL=NULL,start=NULL,maxit=100,trace=0,epsilon1=1e-6,epsilon2=1e-6,zeroN=NULL, digits=max(3,getOption("digits")-3)) {
	if(class(obj)!="readCatdata" && class(obj)!="satMarML") {
		stop("obj is not of the class readCatdata or satMarML")
	}
	if(class(obj)=="readCatdata" && obj$miss==TRUE) {
		stop("With missing data, obj must be of class satMarML")
	}
	tole<-1e-323 #tol argument for solve function

	R<-obj$R
	S<-obj$S
	miss<-obj$miss
	Rp<-obj$Rp
	l<-obj$l
	Tt<-obj$Tt
	Zs<-obj$Zs
	Zbs<-obj$Zbs
	Ns<-obj$Ns
	Nst<-obj$Nst
	Nsm<-obj$Nsm
	nstm<-obj$nstm
	nsmm<-obj$nsmm
	ps<-obj$ps
	pst<-obj$pst
	pbst<-obj$pbst
	bs<-obj$bs
	Bs<-obj$Bs
	b<-obj$b
	B<-obj$B
	theta<-obj$theta
	Vtheta<-obj$Vtheta
	thetab<-obj$thetab
	Vthetab<-obj$Vthetab

	if(miss==TRUE) {
		Zst<-obj$Zst
		Zbst<-obj$Zbst
		missing<-obj$missing
		QvMCAR<-obj$QvMCAR
		alphast<-obj$alphast
	} else {
		missing<-"MCAR"
		alphast<-Nst
		for(s in 1:S) {
			alphast[[s]]<-1
		}
	}

	restr<-kronecker(diag(S),rep(1,R))
	if(is.null(A)) {
		A<-kronecker(diag(S),cbind(diag(R-1),rep(-1,R-1)))
		u<-S*(R-1)
	} else {
		if(!is.matrix(A)) {
			A<-t(A)
		}
		if(ncol(A)!=S*R) {
			stop(paste("The number of columns of A must be equal to S*R=",S*R,".",sep=""))
		}
		u<-nrow(A)
		if(u>S*(R-1)) {
			stop(paste("The number of rows of A must be less than or equal to S*(R-1)=",S*(R-1),".",sep=""))
		}
		if(qr(A)$rank!=u) {
			stop("Rank of A and its number of rows must be the same.")
		}
		if(!all(A%*%restr==matrix(0,u,S))) {
			stop("The lines of A are not orthogonal to the columns of the\nnatural restrictions (kronecker(diag(S),rep(1,R))).")
		}
	}
	if(is.null(X) && is.null(U) && is.null(XL) && is.null(UL)) {
		stop("Specify X or XL for freedom equation formulation or U or UL for constraint formulation.")
	}
	if(  ((!is.null(X)) && (!is.null(U)))   |  ((!is.null(X))  && (!is.null(XL)))  |
	   ((!is.null(X)) && (!is.null(UL)))  |  ((!is.null(U))  && (!is.null(XL)))  |
	   ((!is.null(U)) && (!is.null(UL)))  |  ((!is.null(XL)) && (!is.null(UL)))  ) {
		stop("Specify just X or U or XL or UL, not more than one matrix.")
	}

	#Obtém a matriz U (ou X) a partir da X (ou U), via decomposicao valor singular de uma matriz Q
	#ortogonal a X (ou U). A matriz Q, ortogonal a X, e da forma: Q = Ic - X*inv(X'X)*X'
	#Q possui posto c-r, onde r e o posto de X.
	COMPLEMENT<-function(A) {
		K<-nrow(A)
		J<-ncol(A)
		if(J!=K) { #A nao e' quadrada de posto completo */
			if(J > K) {
				A<-t(A)
				K<-nrow(A)
				J<-ncol(A)
			}
			Q<-diag(K)-A%*%solve(t(A)%*%A,tol=tole)%*%t(A)
			Q.svd<-svd(Q) #decomposicao do valor singular, Q=Q.svd$u %*% diag(Q.svd$d) %*% t(Q.svd$v) ou Q=t(Q.svd$v %*% (t(Q.svd$u) * Q.svd$d))
			#base ortogonal
			B1<-Q.svd$d
			B2<-Q.svd$u
			B3<-Q.svd$v
			CAUX<-matrix(0, nrow=nrow(B2), ncol=1)
			for(i in 1:length(B1)) {
				if(round(B1[i],10)==1) {
					CAUX<-cbind(CAUX,B2[,i])
				}
			}
			if(J < K) {
				CAUX<-CAUX[,2:ncol(CAUX)]
			}
		} else {
			CAUX<-0
		}
		return(CAUX)
	}

	if(is.null(U) && is.null(UL)) {

		if(is.null(XL)) {
			if(!is.matrix(X)) {
				X<-as.matrix(X)
			}
			if(nrow(X)!=S*R) {
				stop(paste("X must have S*R=",S*R," rows.",sep=""))
			}
			p<-ncol(X)
			if(qr(X)$rank!=p) {
				stop("Rank of X and its number of columns must be the same.")
			}
			if(qr(cbind(restr,X))$rank!=(S+p)) {
				stop("The columns of X are not linearly independent from the columns\nof the natural restrictions (kronecker(diag(S),rep(1,R))).")
			}
			XL<-A%*%X
			if(!is.matrix(XL)) {
				XL<-as.matrix(XL)
			}
		} else {
			if(!is.matrix(XL)) {
				XL<-as.matrix(XL)
			}
			if(nrow(XL)!=u) {
				stop("XL and A must have the same number of rows.")
			}
			p<-ncol(XL)
			if(p>u) {
				stop("XL can't have more columns than rows.")
			}
			if(qr(XL)$rank!=p) {
				stop("Rank of XL and its number of columns must be the same.")
			}
		}
		UL<-t(COMPLEMENT(XL))
		form<-"freedom equation"

	} else {

		if(is.null(UL)) {
			if(!is.matrix(U)) {
				U<-t(U)
			}
			if(ncol(U)!=S*R) {
				stop("The number of columns of U have to be the same of the number of S*R.")
			}
			p<-S*(R-1)-nrow(U)
			if(p==0) {
				stop("The number of rows and columns of U must not be the same.")
			}
			if(qr(U)$rank!=(S*(R-1)-p)) {
				stop("Rank of U and its number of rows must be the same.")
			}
			XL<-A%*%COMPLEMENT(rbind(U,t(restr)))
			UL<-t(COMPLEMENT(XL))
		} else {
			if(!is.matrix(UL)) {
				UL<-t(UL)
			}
			if(ncol(UL)!=u) {
				stop("The number of columns of UL have to be the same of the number of rows of A.")
			}
			p<-u-nrow(UL)
			if(p==0) {
				stop("The number of rows and columns of UL must not be the same.")
			}
			if(qr(UL)$rank!=(u-p)) {
				stop("Rank of UL and its number of rows must be the same.")
			}
			XL<-COMPLEMENT(UL)
			if(!is.matrix(XL)) {
				XL<-as.matrix(XL)
			}
		}
		form<-"constraint"

	}
	if(all(XL==0)) {
		X<-0
	} else {
		X<-t(A)%*%solve(A%*%t(A),tol=tole)%*%XL
	}
	Aa<-A
	Xa<-X
	XLa<-XL
	if(u<S*(R-1)) {
		A0<-t(COMPLEMENT(rbind(A,t(restr))))
		Aa<-rbind(A,A0)
		if(all(XL==0)) {
			Xa<-cbind(t(A0)%*%solve(A0%*%t(A0),tol=tole))
			XLa<-rbind( matrix(0,u,S*(R-1)-u), diag(S*(R-1)-u) )
		} else {
			Xa<-cbind(X,t(A0)%*%solve(A0%*%t(A0),tol=tole))
			XLa<-rbind( cbind(XL,matrix(0,u,S*(R-1)-u)), cbind(matrix(0,S*(R-1)-u,p),diag(S*(R-1)-u)) )
		}
	}
	if(!is.null(zeroN)) {
		if(!is.matrix(zeroN) || (nrow(zeroN)!=S) || (ncol(zeroN)!=max(Tt))) {
			if(is.vector(zeroN,mode="numeric")) {
				zeroN<-matrix(zeroN[1],S,max(Tt))
			} else {
				stop("zeroN must be a number or a matrix with dimensions S x (max(Tt))")
			}
		}
		if(any(zeroN<0 | zeroN>0.5)) {
			stop("The elements of zeroN must be nonnegative and less than or equal 0.5.")
		}
	} else {
		if(miss==TRUE) {
			zeroN<-1/(cbind(rep(R,S),Rp)*nstm)
			zeroN[zeroN==Inf]<-0
		} else {
			zeroN<-1/(rep(R,S)*nstm)
		}
	}
	if(!is.numeric(maxit) || maxit<1) {
		stop("maxit must be >= 1")
	}
	if(!any(trace==c(0,1,2))) {
		stop("trace must be 0, 1 or 2")
	}
	if(!is.numeric(epsilon1) || epsilon1<=0) {
		stop("epsilon1 must be > 0")
	}
	if(!is.numeric(epsilon2) || epsilon2<=0 || epsilon2>=0.1) {
		stop("epsilon2 must be > 0 and < 0.1")
	}

	#Score vector S_{1LL}(beta)
	S1LL<-function(thetaH,R,S,Tt,Zst,Nst,nsmm,p,u,Xa) {
		s1ll<-numeric(p+S*(R-1)-u)
		for(s in 1:S) {
			thetaHs<-thetaH[((s-1)*R+1):(s*R)]
			aux<-Nst[[paste("st",s,".",1,sep="")]]-nsmm[s]*thetaHs
			if(Tt[s]>1) {
				for(tt in 2:Tt[s]) {
					zst<-Zst[[paste("st",s,".",tt,sep="")]]
					aux<-aux+c( diag(thetaHs)%*%zst%*%solve(diag(c(t(zst)%*%thetaHs)),tol=tole)%*%Nst[[paste("st",s,".",tt,sep="")]] )
				}
			}
			s1ll<-s1ll+c(t(Xa[((s-1)*R+1):(s*R),])%*%aux)
		}
		return(s1ll)
	}

	#Fisher information matrix i_{1LL}(beta,\{alpha_{st}^M\})
	IF1LL<-function(thetaH,missing,R,S,Tt,Rp,Zst,nsmm,alphast,p,u,Xa) {
		if1ll<-matrix(0,p+S*(R-1)-u,p+S*(R-1)-u)
		for(s in 1:S) {
			thetaHs<-thetaH[((s-1)*R+1):(s*R)]
			aux<-diag(R)
			if(Tt[s]>1) {
				for(tt in 2:Tt[s]) {
					zst<-Zst[[paste("st",s,".",tt,sep="")]]
					if(missing=="MAR") {
						alpha<-alphast[[paste("st",s,".",tt,sep="")]]
						aux<-aux-diag(c(zst%*%alpha))+diag(c(diag(thetaHs)%*%zst%*%solve(diag(c(t(zst)%*%thetaHs)),tol=tole)%*%alpha))%*%zst%*%t(zst)
					}
					if(missing=="MCAR") {
						alpha<-alphast[[paste("st",s,".",tt,sep="")]]
						aux<-aux-diag(rep(alpha,R))+diag(c(diag(thetaHs)%*%zst%*%solve(diag(c(t(zst)%*%thetaHs)),tol=tole)%*%rep(alpha,Rp[s,tt-1])))%*%zst%*%t(zst)
					}
				}
			}
			if1ll<-if1ll+nsmm[s]*t(Xa[((s-1)*R+1):(s*R),])%*%aux%*%(diag(thetaHs)-thetaHs%*%t(thetaHs))%*%Xa[((s-1)*R+1):(s*R),]
		}
		return(if1ll)
	}

	#Likelihood ratio statistic of H given MAR or any missing mechanism contained in MAR
	Qv<-function(thetaH,theta,R,S,Zs,Ns) {
		Qv<-0
		for(s in 1:S) {
			Qv<-Qv-2*c( t(Ns[[paste("s",s,sep="")]][Ns[[paste("s",s,sep="")]]!=0]) %*%
				   (log(c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)])[Ns[[paste("s",s,sep="")]]!=0])
				    -log(c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)])[Ns[[paste("s",s,sep="")]]!=0]) ) )
		}
		return(Qv)
	}
	conv<-FALSE
	it<-0
	theta0<-theta
	Vtheta0<-Vtheta
	if(miss==FALSE) {
		for(s in 1:S) {
			if(any(theta[((s-1)*R+1):(s*R)]==0)) {
				theta0[((s-1)*R+1):(s*R)]<-Nst[[paste("st",s,".",1,sep="")]]
				theta0[((s-1)*R+1):(s*R)][theta0[((s-1)*R+1):(s*R)]==0]<-zeroN[s,1]
				theta0[((s-1)*R+1):(s*R)]<-theta0[((s-1)*R+1):(s*R)]/sum(theta0[((s-1)*R+1):(s*R)])
				Vtheta0[((s-1)*R+1):(s*R),((s-1)*R+1):(s*R)]<-(diag(theta0[((s-1)*R+1):(s*R)])-theta0[((s-1)*R+1):(s*R)]%*%t(theta0[((s-1)*R+1):(s*R)]))/nsmm[s]
			}
		}
	}
	if(is.null(start)) {
		aux<-solve(Aa%*%solve(diag(theta0))%*%Vtheta0%*%solve(diag(theta0))%*%t(Aa),tol=tole)
		beta<-c( solve(t(XLa)%*%aux%*%XLa,tol=tole)%*%t(XLa)%*%aux%*%Aa%*%log(theta0) )
	} else {
		if(!is.vector(start,mode="numeric")) {
			stop("start must be a numeric vector")
		}
		if(u==S*(R-1)) {
			if(form=="constraint") {
				stop("start may be used just with the freedom equation formulation")
			}
			if(length(start)!=p) {
				stop("start must have the same length of the number of parameters")
			}
			beta<-start
		} else {
			stop("start may be supplied just if the number of rows of A is S(R-1)")
		}
	}
	thetaH<-c(solve(diag(c(kronecker(diag(S),matrix(1,R,R))%*%exp(Xa%*%beta))),tol=tole)%*%exp(Xa%*%beta))
	QvH<-Qv(thetaH,theta,R,S,Zs,Ns)
	if(trace==1) {
		cat("Iteration - Likelihood ratio statistic of the log-linear model\n")
		cat(it," - ",QvH,"\n")
	}
	if(trace==2) {
		cat("Iteration ",it," - LRS(loglin) ",QvH," - Parameter estimates:\n",sep="")
		print.default(format(round(beta,digits=digits),digits=digits),quote=FALSE,print.gap=2)
		cat("\n")
	}
	while((conv==FALSE) & (it<maxit)) {
		it<-it+1
		betaa<-beta
		QvHa<-QvH
		beta<-betaa+c( solve(IF1LL(thetaH,missing,R,S,Tt,Rp,Zst,nsmm,alphast,p,u,Xa),tol=tole) %*% S1LL(thetaH,R,S,Tt,Zst,Nst,nsmm,p,u,Xa) )
		thetaH<-c(solve(diag(c(kronecker(diag(S),matrix(1,R,R))%*%exp(Xa%*%beta))),tol=tole)%*%exp(Xa%*%beta))
		if(any((thetaH<=0) | (thetaH>=1))) {
			if(any((thetaH<0) | (thetaH>1))) {
				stop(paste("Any of the estimated probabilities obtained by the iterative\nprocess are outside the parameter space. (iteration ",it,")",sep=""))
			} else {
				if(min(thetaH)==0) {
					warning(paste("Any of the estimated probabilities obtained by the iterative process are in the\nborderline of the parameter space (are equal 0). (iteration ",it,")",sep=""))
				}
				if(max(thetaH)==1) {
					warning(paste("Any of the estimated probabilities obtained by the iterative process are in the\nborderline of the parameter space (are equal 1). (iteration ",it,")",sep=""))
				}
			}
		}
		QvH<-Qv(thetaH,theta,R,S,Zs,Ns)
		if(trace==1) {
			cat(it," - ",QvH,"\n")
		}
		if(trace==2) {
			cat("Iteration ",it," - LRS(loglin) ",QvH," - Parameter estimates:\n",sep="")
			print.default(format(round(beta,digits=digits),digits=digits),quote=FALSE,print.gap=2)
			cat("\n")
		}
		if(abs(QvHa-QvH)<epsilon1 && max(abs(beta-betaa))<epsilon2) {
			conv<-TRUE
		}
	}
	if(conv==FALSE) {
		warning(paste("\nThe iterative process has not attained the convergence criterions with ",maxit," iterations.\n",sep=""))
	}
	Vbeta<-solve(IF1LL(thetaH,missing,R,S,Tt,Rp,Zst,nsmm,alphast,p,u,Xa),tol=tole)
	VLL<-matrix(0,S*R,S*R)
	for(s in 1:S) {
		VLL[((s-1)*R+1):(s*R),((s-1)*R+1):(s*R)]<-diag(thetaH[((s-1)*R+1):(s*R)])-thetaH[((s-1)*R+1):(s*R)]%*%t(thetaH[((s-1)*R+1):(s*R)])
	}
	VthetaH<-VLL%*%Xa%*%Vbeta%*%t(Xa)%*%VLL
	thetabH<-c(kronecker(diag(S),cbind(diag(R-1),rep(0,R-1)))%*%thetaH)
	VthetabH<-kronecker(diag(S),cbind(diag(R-1),rep(0,R-1)))%*%VthetaH%*%t(kronecker(diag(S),cbind(diag(R-1),rep(0,R-1))))
	if(u<S*(R-1)) {
		if(p>0) {
			beta<-beta[1:p]
			Vbeta<-Vbeta[1:p,1:p]
		} else {
			beta<-0
			Vbeta<-0
		}
	}
	Fu<-c(A%*%log(theta0)) #F=A*log(theta) unrestricted
	FH<-c(XL%*%beta)  #F=A*log(theta) under log-linear model
	VFu<-A%*%solve(diag(theta0),tol=tole)%*%Vtheta0%*%solve(diag(theta0),tol=tole)%*%t(A)
	VFH<-XL%*%Vbeta%*%t(XL)
	if(miss==TRUE) {
		QvHMCAR<-QvH+QvMCAR #Likelihood ratio statistic of H and MCAR given MAR
		QpHMCAR<-0 #Pearson statistic of H and MCAR given MAR
		QnHMCAR<-0 #Neyman statistic of H and MCAR given MAR
	}
	QpH<-0 #Pearson statistic of H given missing mechanism (MAR or MCAR)
	QnH<-0 #Neyman statistic of H given missing mechanism (MAR or MCAR)
	for(s in 1:S) {
		if(miss==TRUE) {
			pp<-pst[[paste("st",s,".",1,sep="")]]
			if(any(pp==0)) {
				pp<-Nst[[paste("st",s,".",1,sep="")]]
				pp[pp==0]<-zeroN[s,1]
				pp<-pp/sum(pp)
			}
			if(missing!="MAR") {
				th<-c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)])
			}
		} else {
			th<-pst[[paste("st",s,".",1,sep="")]]
			if(any(th==0)) {
				th<-Nst[[paste("st",s,".",1,sep="")]]
				th[th==0]<-zeroN[s,1]
				th<-th/sum(th)
			}
		}
		if(Tt[s]>1) {
			for(tt in 2:Tt[s]) {
				pp1<-pst[[paste("st",s,".",tt,sep="")]]
				if(any(pp1==0)) {
					pp1<-Nst[[paste("st",s,".",tt,sep="")]]
					pp1[pp1==0]<-zeroN[s,tt]
					pp1<-pp1/sum(pp1)
				}
				pp<-c(pp,pp1)
			}
		}
		if(missing=="MAR") {
			QpH<-QpH+c( t(c(t(Zs[[paste("s",s,sep="")]])%*%(theta[((s-1)*R+1):(s*R)]-thetaH[((s-1)*R+1):(s*R)])))%*%
				   (diag(Ns[[paste("s",s,sep="")]])%*%solve(diag(c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)])),tol=tole)%*%solve(diag(c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)])),tol=tole))%*%
			(c(t(Zs[[paste("s",s,sep="")]])%*%(theta[((s-1)*R+1):(s*R)]-thetaH[((s-1)*R+1):(s*R)]))) )
			QnH<-QnH+c( t(rep(1,R+l[s])-c(solve(diag(c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)])),tol=tole)%*%t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)]))%*%
				   diag(Ns[[paste("s",s,sep="")]])%*%
			(rep(1,R+l[s])-c(solve(diag(c(t(Zs[[paste("s",s,sep="")]])%*%theta[((s-1)*R+1):(s*R)])),tol=tole)%*%t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)])) )
		} else {
			QpH<-QpH+c( t(c(t(Zs[[paste("s",s,sep="")]])%*%(theta[((s-1)*R+1):(s*R)]-thetaH[((s-1)*R+1):(s*R)])))%*%
				   (diag(Nsm[[paste("s",s,sep="")]])%*%solve(diag(c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)])),tol=tole))%*%
			(c(t(Zs[[paste("s",s,sep="")]])%*%(theta[((s-1)*R+1):(s*R)]-thetaH[((s-1)*R+1):(s*R)]))) )
			QnH<-QnH+c( t(c(t(Zs[[paste("s",s,sep="")]])%*%(theta[((s-1)*R+1):(s*R)]-thetaH[((s-1)*R+1):(s*R)])))%*%
				   (diag(Nsm[[paste("s",s,sep="")]])%*%solve(diag(th),tol=tole))%*%
			(c(t(Zs[[paste("s",s,sep="")]])%*%(theta[((s-1)*R+1):(s*R)]-thetaH[((s-1)*R+1):(s*R)]))) )
		}
		if(miss==TRUE) {
			QpHMCAR<-QpHMCAR+c( t(ps[[paste("s",s,sep="")]]-c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)]))%*%
					   (diag(Nsm[[paste("s",s,sep="")]])%*%solve(diag(c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)])),tol=tole))%*%
			(ps[[paste("s",s,sep="")]]-c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)])) )
			QnHMCAR<-QnHMCAR+c( t(ps[[paste("s",s,sep="")]]-c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)]))%*%
					   (diag(Nsm[[paste("s",s,sep="")]])%*%solve(diag(pp),tol=tole))%*%
			(ps[[paste("s",s,sep="")]]-c(t(Zs[[paste("s",s,sep="")]])%*%thetaH[((s-1)*R+1):(s*R)])) )
		}
	}
	glH<-u-p
	glHMCAR<-glH+S+sum(l)-sum(Tt)
	QwH<-0
	if(glH>0) {
		QwH<-c(t(UL%*%Fu)%*%solve(UL%*%VFu%*%t(UL),tol=tole)%*%UL%*%Fu)
	}
	ystH<-Nst
	for(s in 1:S) {
		for(tt in 1:Tt[s]) {
			if(missing=="MAR") {
				if(tt==1) {
					zst<-diag(R)
				} else {
					zst<-Zst[[paste("st",s,".",tt,sep="")]]
				}
				ystH[[paste("st",s,".",tt,sep="")]]<-c( nsmm[s]*diag(thetaH[((s-1)*R+1):(s*R)])%*%zst%*%alphast[[paste("st",s,".",tt,sep="")]] )
			}
			if(missing=="MCAR") {
				ystH[[paste("st",s,".",tt,sep="")]]<-nstm[s,tt]*thetaH[((s-1)*R+1):(s*R)]
			}
		}
	}
	if(miss==TRUE) {
		res<-list(obj=obj,missing=missing,R=R,S=S,miss=miss,Rp=Rp,l=l,Tt=Tt,
			  Zs=Zs,Zst=Zst,Zbs=Zbs,Zbst=Zbst,Ns=Ns,Nst=Nst,Nsm=Nsm,nstm=nstm,nsmm=nsmm,ps=ps,pst=pst,pbst=pbst,bs=bs,Bs=Bs,b=b,B=B,
			  theta=theta,Vtheta=Vtheta,thetab=thetab,Vthetab=Vthetab,alphast=alphast,maxit=maxit,epsilon1=epsilon1,
			  epsilon2=epsilon2,zeroN=zeroN,call=match.call(),A=A,X=X,U=U,XL=XL,UL=UL,Aa=Aa,Xa=Xa,XLa=XLa,form=form,p=p,u=u,conv=conv,
			  it=it,thetaH=thetaH,VthetaH=VthetaH,thetabH=thetabH,VthetabH=VthetabH,beta=beta,Vbeta=Vbeta,Fu=Fu,VFu=VFu,FH=FH,VFH=VFH,
			  QvH=QvH,QpH=QpH,QnH=QnH,QwH=QwH,glH=glH,QvHMCAR=QvHMCAR,QpHMCAR=QpHMCAR,QnHMCAR=QnHMCAR,glHMCAR=glHMCAR,ystH=ystH)
	} else {
		res<-list(obj=obj,R=R,S=S,miss=miss,Rp=Rp,l=l,Tt=Tt,
			  Zs=Zs,Zbs=Zbs,Ns=Ns,Nst=Nst,Nsm=Nsm,nstm=nstm,nsmm=nsmm,ps=ps,pst=pst,pbst=pbst,bs=bs,Bs=Bs,b=b,B=B,
			  theta=theta,Vtheta=Vtheta,thetab=thetab,Vthetab=Vthetab,maxit=maxit,epsilon1=epsilon1,
			  epsilon2=epsilon2,zeroN=zeroN,call=match.call(),A=A,X=X,U=U,XL=XL,UL=UL,Aa=Aa,Xa=Xa,XLa=XLa,form=form,p=p,u=u,conv=conv,
			  it=it,thetaH=thetaH,VthetaH=VthetaH,thetabH=thetabH,VthetabH=VthetabH,beta=beta,Vbeta=Vbeta,Fu=Fu,VFu=VFu,FH=FH,VFH=VFH,
			  QvH=QvH,QpH=QpH,QnH=QnH,QwH=QwH,glH=glH,ystH=ystH)
	}
	class(res)<-"loglinML"
	res
}


print.loglinML<-function(x, digits=max(3,getOption("digits")-3), ...) {
	cat("\nCall: ",deparse(x$call),"\n",sep="")

	if(x$form=="freedom equation") {
		pp<-matrix(0,x$p,4)
		pp[,1]<-x$beta
		pp[,2]<-sqrt(diag(as.matrix(x$Vbeta)))
		pp[,3]<-pp[,1]/pp[,2]
		pp[,4]<-1-pchisq((pp[,3])^2,1)
		colnames(pp)<-c("estimate","std.error","z-value","p-value")
		cat("\nMaximum likelihood estimates of the parameters of the log-linear model",sep="")
		if(x$miss==TRUE) {
			cat(" under ",x$missing,sep="")
		}
		cat(":\n",sep="")
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)
	}

	pp<-matrix(0,4,2)
	pp[1,]<-c(x$QvH,ifelse(x$glH>0,1-pchisq(x$QvH,x$glH),1))
	pp[2,]<-c(x$QpH,ifelse(x$glH>0,1-pchisq(x$QpH,x$glH),1))
	pp[3,]<-c(x$QnH,ifelse(x$glH>0,1-pchisq(x$QnH,x$glH),1))
	pp[4,]<-c(x$QwH,ifelse(x$glH>0,1-pchisq(x$QwH,x$glH),1))
	rownames(pp)<-c("Likelihood ratio","Pearson","Neyman","Wald")
	colnames(pp)<-c("statistic","p-value")

	cat("\nGoodness of fit of the log-linear model",sep="")
	if(x$miss==TRUE) {
		cat(" given ",x$missing,sep="")
	}
	cat(" (d.f.=",x$glH,"):\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)

	if(x$miss==TRUE) {
		pp<-matrix(0,3,2)
		pp[1,]<-c(x$QvHMCAR,1-pchisq(x$QvHMCAR,x$glHMCAR))
		pp[2,]<-c(x$QpHMCAR,1-pchisq(x$QpHMCAR,x$glHMCAR))
		pp[3,]<-c(x$QnHMCAR,1-pchisq(x$QnHMCAR,x$glHMCAR))
		rownames(pp)<-c("Likelihood ratio","Pearson","Neyman")
		colnames(pp)<-c("statistic","p-value")
		cat("\nGoodness of fit of the log-linear model and MCAR given MAR (d.f.=",x$glHMCAR,"):\n",sep="")
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)
	}

	cat("\n")
}

summary.loglinML<-function(object, digits=max(3,getOption("digits")-3), ...) {

	cat("\nCall: ",deparse(object$call),"\n\n",sep="")
	pp<-matrix(0,object$S,object$R)
	epp<-matrix(0,object$S,object$R)
	for(s in 1:object$S) {
		pp[s,]<-object$thetaH[((s-1)*object$R+1):(s*object$R)]
		epp[s,]<-sqrt(diag(object$VthetaH[((s-1)*object$R+1):(s*object$R),((s-1)*object$R+1):(s*object$R)]))
	}
	cat("\nMaximum likelihood estimates of the probabilities under the log-linear model (LLM):\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	cat("\nStandard errors",sep="")
	if(object$miss==TRUE) {
		cat(" (",object$missing,")",sep="")
	}
	cat(":\n",sep="")
	print.default(format(round(epp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	pp<-matrix(0,object$u,4)
	pp[,1]<-object$Fu
	pp[,2]<-sqrt(diag(as.matrix(object$VFu)))
	pp[,3]<-object$FH
	pp[,4]<-sqrt(diag(as.matrix(object$VFH)))
	colnames(pp)<-c("observed","std.error","under the LLM","std.error")
	cat("\n\nMaximum likelihood estimates of the log-linear functions:\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	if(object$form=="freedom equation") {
		pp<-matrix(0,object$p,4)
		pp[,1]<-object$beta
		pp[,2]<-sqrt(diag(as.matrix(object$Vbeta)))
		pp[,3]<-pp[,1]/pp[,2]
		pp[,4]<-1-pchisq((pp[,3])^2,1)
		colnames(pp)<-c("estimate","std.error","z-value","p-value")
		cat("\n\nMaximum likelihood estimates of the parameters of the log-linear model",sep="")
		if(object$miss==TRUE) {
			cat(" under ",object$missing,sep="")
		}
		cat(":\n",sep="")
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)
	}

	if(object$conv==FALSE) {
		cat("\n\nFisher scoring has NOT attained the convergence criterion with ",object$maxit," iterations.\n",sep="")
	} else {
		cat("\n\nFisher scoring attained the convergence criterion in ",object$it," iterations.\n",sep="")
	}

	pp<-matrix(0,4,2)
	pp[1,]<-c(object$QvH,ifelse(object$glH>0,1-pchisq(object$QvH,object$glH),1))
	pp[2,]<-c(object$QpH,ifelse(object$glH>0,1-pchisq(object$QpH,object$glH),1))
	pp[3,]<-c(object$QnH,ifelse(object$glH>0,1-pchisq(object$QnH,object$glH),1))
	pp[4,]<-c(object$QwH,ifelse(object$glH>0,1-pchisq(object$QwH,object$glH),1))
	rownames(pp)<-c("Likelihood ratio","Pearson","Neyman","Wald")
	colnames(pp)<-c("statistic","p-value")

	cat("\n\nGoodness of fit of the log-linear model",sep="")
	if(object$miss==TRUE) {
		cat(" given ",object$missing,sep="")
	}
	cat(" (d.f.=",object$glH,"):\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)

	if(object$miss==TRUE) {
		pp<-matrix(0,3,2)
		pp[1,]<-c(object$QvHMCAR,1-pchisq(object$QvHMCAR,object$glHMCAR))
		pp[2,]<-c(object$QpHMCAR,1-pchisq(object$QpHMCAR,object$glHMCAR))
		pp[3,]<-c(object$QnHMCAR,1-pchisq(object$QnHMCAR,object$glHMCAR))
		rownames(pp)<-c("Likelihood ratio","Pearson","Neyman")
		colnames(pp)<-c("statistic","p-value")
		cat("\n\nGoodness of fit of the log-linear model and MCAR given MAR (d.f.=",object$glHMCAR,"):\n",sep="")
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)
	}

	if(object$miss==TRUE) {
		cat("\n\nAugmented estimated frequencies under log-linear model and ",object$missing,":\n",sep="")
		for(s in 1:object$S) {
			if(object$S>1) {
				cat("\nSubpopulation ",s,"\n",sep="")
			}
			pp<-matrix(0,object$Tt[s],object$R)
			pp[1,]<-object$yst[[paste("st",s,".",1,sep="")]]
			if(object$Tt[s]>1) {
				for(tt in 2:object$Tt[s]) {
					pp[tt,]<-object$yst[[paste("st",s,".",tt,sep="")]]
				}
			}
			print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)
		}
	} else {
		cat("\n\nEstimated frequencies under log-linear model:\n",sep="")
		pp<-matrix(0,object$S,object$R)
		for(s in 1:object$S) {
			pp[s,]<-object$yst[[paste("st",s,".",1,sep="")]]
		}
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)
	}

	cat("\n")
}





#Fit functional linear models by weighted least squares to the marginal probabilities of categorization
funlinWLS<-function(model,obj=NULL,theta=NULL,Vtheta=NULL,A1=NULL,A2=NULL,A3=NULL,A4=NULL,A5=NULL,
																				A6=NULL,A7=NULL,A8=NULL,A9=NULL,X=NULL,U=NULL,XL=NULL,UL=NULL,zeroN=NULL,PI1=NULL,
																				PI2=NULL,PI3=NULL,PI4=NULL,PI5=NULL,PI6=NULL,PI7=NULL,PI8=NULL,PI9=NULL) {
	if(!is.vector(model,mode="character")) {
		stop("model must be a character vector")
	}
	if(any(model!="lin" & model!="log" & model!="exp" & model!="add")) {
		stop("model must have only elements equal 'lin', 'log', 'exp' or 'add'")
	}
	tole<-1e-323 #tol argument for solve function

	nfunc<-length(model) #Total number of functions

	if( (!is.null(obj)) && (  (!is.null(theta)) | (!is.null(Vtheta))  ) ) {
		stop("Supply just obj or theta and Vtheta")
	}
	if(is.null(obj)) {
		if( is.null(theta) | is.null(Vtheta)  ) {
			stop("obj or theta and Vtheta must be supplied")
		}
		if(!is.vector(theta,mode="numeric")) {
			stop("theta must be a numeric vector")
		}
		if(!is.matrix(Vtheta)) {
			stop("Vtheta must be a matrix")
		}
		if(ncol(Vtheta)!=nrow(Vtheta)) {
			stop("Vtheta must be a square matrix")
		}
		if(ncol(Vtheta)!=length(theta)) {
			stop("Vtheta must have the same dimension of the length of theta")
		}
		Fu<-theta #F unrestricted
		VFu<-Vtheta
		type<-"functional linear (out)"
	} else {
		if(!any(class(obj)==c("readCatdata","satMarML","satMcarWLS"))) {
			stop("obj is not of the class readCatdata, satMarML or satMcarWLS")
		}
		if(class(obj)=="readCatdata" && obj$miss==TRUE) {
			stop("With missing data, obj must be of class satMarML or satMcarWLS")
		}
		theta<-obj$theta #F unrestricted
		Vtheta<-obj$Vtheta
		Fu<-theta #F unrestricted
		VFu<-Vtheta
		R<-obj$R
		S<-obj$S
		miss<-obj$miss
		Tt<-obj$Tt
		b<-obj$b
		B<-obj$B
		restr<-kronecker(diag(S),rep(1,R))
		type<-"functional linear (obj)"
		if(nfunc==1 && model=="lin") {
			type<-"linear"
			if(is.null(A1)) {
				A1<-kronecker(diag(S),cbind(diag(R-1),rep(0,R-1)))
				u<-S*(R-1)
			} else {
				if(!is.matrix(A1)) {
					A1<-t(A1)
				}
				if(ncol(A1)!=S*R) {
					stop(paste("The number of columns of A1 must be equal to S*R=",S*R,".",sep=""))
				}
				u<-nrow(A1)
				if(u>S*(R-1)) {
					stop(paste("The number of rows of A1 must be less than or equal to S*(R-1)=",S*(R-1),".",sep=""))
				}
				if(qr(A1)$rank!=u) {
					stop("Rank of A1 and its number of rows must be the same.")
				}
				if(qr(cbind(t(A1),restr))$rank!=(u+S)) {
					stop("The lines of A1 are not linearly independent from the columns\nof the natural restrictions (kronecker(diag(S),rep(1,R))).")
				}
			}
		}
		if(nfunc==2 && all(model==c("lin","log"))) {
			if(is.null(A1)) {
				A1<-kronecker(diag(S),cbind(diag(R-1),rep(-1,R-1)))
				u<-S*(R-1)
				type<-"log-linear"
			} else {
				if(!is.matrix(A1)) {
					A1<-t(A1)
				}
				if(ncol(A1)!=S*R) {
					stop(paste("The number of columns of A1 must be equal to S*R=",S*R,".",sep=""))
				}
				u<-nrow(A1)
				if(u>S*(R-1)) {
					stop(paste("The number of rows of A1 must be less than or equal to S*(R-1)=",S*(R-1),".",sep=""))
				}
				if(qr(A1)$rank!=u) {
					stop("Rank of A1 and its number of rows must be the same.")
				}
				if(all(A1%*%restr==matrix(0,u,S))) {
					type<-"log-linear"
				}
			}
		}
	}

	if(type!="functional linear (out)" && miss==FALSE) {
		if(!is.null(zeroN)) {
			if(!is.vector(zeroN,mode="numeric")) {
				stop("zeroN must be a number or a vector with length equal S")
			}
			if(length(zeroN)!=S) {
				zeroN<-rep(zeroN[1],S)
			}
			if(any(zeroN<0 | zeroN>0.5)) {
				stop("The elements of zeroN must be nonnegative and less than or equal 0.5.")
			}
		} else {
			zeroN<-1/(rep(R,S)*obj$nsmm)
		}
		for(s in 1:S) {
			if(any(Fu[((s-1)*R+1):(s*R)]==0)) {
				Fu[((s-1)*R+1):(s*R)]<-obj$Nst[[paste("st",s,".",1,sep="")]]
				Fu[((s-1)*R+1):(s*R)][Fu[((s-1)*R+1):(s*R)]==0]<-zeroN[s]
				Fu[((s-1)*R+1):(s*R)]<-Fu[((s-1)*R+1):(s*R)]/sum(Fu[((s-1)*R+1):(s*R)])
				VFu[((s-1)*R+1):(s*R),((s-1)*R+1):(s*R)]<-(diag(Fu[((s-1)*R+1):(s*R)])-Fu[((s-1)*R+1):(s*R)]%*%t(Fu[((s-1)*R+1):(s*R)]))/obj$nsmm[s]
			}
		}
		theta<-Fu
		Vtheta<-VFu
	}

	iA<-0
	iPI<-0
	for(i in nfunc:1) {
		if(model[i]=="lin") {
			iA<-iA+1
			if(iA>9) {
				stop("You can apply a maximum of 9 linear operations.")
			}
			A<-eval(as.name(paste("A",iA,sep="")))
			if(!is.matrix(A)) {
				stop(paste("A",iA," must be a matrix.",sep=""))
			}
			if(ncol(A)!=length(Fu)) {
				stop(paste("Matrix A",iA," should have ",length(Fu)," columns. Check this matrix and/or previous operations.",sep=""))
			}
			Fu<-c( A %*% Fu )
			VFu<-A %*% VFu %*% t(A)
		}
		if(model[i]=="log") {
			if(any(Fu<=0)) {
				if(i==nfunc) {
					stop("It is not possible to apply log operation because there is at least one theta that is nonpositive.")
				} else {
					stop(paste("It is not possible to apply log operation to the result of the last ",nfunc-i," operation(s)\nspecified in the vector model, from right to left, because there is at least one number that is nonpositive.",sep=""))
				}
			}
			A<-solve(diag(Fu,nrow=length(Fu)),tol=tole)
			Fu<-log( Fu )
			VFu<-A %*% VFu %*% t(A)
		}
		if(model[i]=="exp") {
			Fu<-exp( Fu )
			A<-diag(Fu,nrow=length(Fu))
			VFu<-A %*% VFu %*% t(A)
		}
		if(model[i]=="add") {
			iPI<-iPI+1
			if(iPI>9) {
				stop("You can apply a maximum of 9 sum of constants.")
			}
			PI<-eval(as.name(paste("PI",iPI,sep="")))
			if(!is.vector(PI,mode="numeric")) {
				stop(paste("PI",iPI," must be a numeric vector.",sep=""))
			}
			if(length(PI)!=length(Fu)) {
				stop(paste("Vector PI",iPI," should have length ",length(Fu),". Check this vector and/or previous operations.",sep=""))
			}
			Fu<-c( Fu + PI )
		}
	}
	u<-length(Fu)
	if(qr(VFu)$rank!=u) {
		stop(paste("The rank of the observed covariance matrix of your specified functions is ",qr(VFu)$rank,",\nbut it should be equal the length of the observed functions (=",u,").",sep=""))
	}

	#Obtém a matriz U (ou X) a partir da X (ou U), via decomposicao valor singular de uma matriz Q
	#ortogonal a X (ou U). A matriz Q, ortogonal a X, e da forma: Q = Ic - X*inv(X'X)*X'
	#Q possui posto c-r, onde r e o posto de X.
	COMPLEMENT<-function(A) {
		K<-nrow(A)
		J<-ncol(A)
		if(J!=K) { #A nao e' quadrada de posto completo */
			if(J > K) {
				A<-t(A)
				K<-nrow(A)
				J<-ncol(A)
			}
			Q<-diag(K)-A%*%solve(t(A)%*%A,tol=tole)%*%t(A)
			Q.svd<-svd(Q) #decomposicao do valor singular, Q=Q.svd$u %*% diag(Q.svd$d) %*% t(Q.svd$v) ou Q=t(Q.svd$v %*% (t(Q.svd$u) * Q.svd$d))
			#base ortogonal
			B1<-Q.svd$d
			B2<-Q.svd$u
			B3<-Q.svd$v
			CAUX<-matrix(0, nrow=nrow(B2), ncol=1)
			for(i in 1:length(B1)) {
				if(round(B1[i],10)==1) {
					CAUX<-cbind(CAUX,B2[,i])
				}
			}
			if(J < K) {
				CAUX<-CAUX[,2:ncol(CAUX)]
			}
		} else {
			CAUX<-0
		}
		return(CAUX)
	}

	if(type=="linear") {


		if(is.null(X) && is.null(U)) {
			stop("Specify X for freedom equation formulation or U for constraint formulation.")
		}
		if((!is.null(X)) && (!is.null(U))) {
			stop("Specify just X or U, not both.")
		}
		if(is.null(U)) {
			if(!is.matrix(X)) {
				X<-as.matrix(X)
			}
			if(nrow(X)!=u) {
				stop("X and A1 must have the same number of rows.")
			}
			p<-ncol(X)
			if(p>u) {
				stop("X can't have more columns than rows.")
			}
			if(qr(X)$rank!=p) {
				stop("Rank of X and its number of columns must be the same.")
			}
			U<-t(COMPLEMENT(X))
			form<-"freedom equation"
			if(!is.matrix(U)) {
				U<-t(U)
			}
		} else {
			if(!is.matrix(U)) {
				U<-t(U)
			}
			if(ncol(U)!=u) {
				stop("The number of columns of U have to be the same of the number of rows of A1.")
			}
			p<-u-nrow(U)
			if(p==0) {
				stop("The number of rows and columns of U must not be the same.")
			}
			if(qr(U)$rank!=(u-p)) {
				stop("Rank of U and its number of rows must be the same.")
			}
			X<-COMPLEMENT(U)
			form<-"constraint"
			if(!is.matrix(X)) {
				X<-as.matrix(X)
			}
		}
		if(u<S*(R-1)) {
			A0<-t(COMPLEMENT(rbind(A1,t(restr))))
			Wpart<-kronecker(diag(S),cbind(diag(R-1),rep(0,R-1)))%*%solve(rbind(A1,t(restr),A0),tol=tole)
			W<-Wpart%*%rbind(cbind(X,matrix(0,u,S*(R-1)-u)),matrix(0,S,S*(R-1)-u+p),cbind(matrix(0,S*(R-1)-u,p),diag(S*(R-1)-u)))
		} else {
			Wpart<-kronecker(diag(S),cbind(diag(R-1),rep(0,R-1)))%*%solve(rbind(A1,t(restr)),tol=tole)
			W<-Wpart%*%rbind(X,matrix(0,S,p))
		}


	}

	if(type=="log-linear") {


		if(is.null(X) && is.null(U) && is.null(XL) && is.null(UL)) {
			stop("Specify X or XL for freedom equation formulation or U or UL for constraint formulation.")
		}
		if(  ((!is.null(X)) && (!is.null(U)))   |  ((!is.null(X))  && (!is.null(XL)))  |
		   ((!is.null(X)) && (!is.null(UL)))  |  ((!is.null(U))  && (!is.null(XL)))  |
		   ((!is.null(U)) && (!is.null(UL)))  |  ((!is.null(XL)) && (!is.null(UL)))  ) {
			stop("Specify just X or U or XL or UL, not more than one matrix.")
		}

		if(is.null(U) && is.null(UL)) {

			if(is.null(XL)) {
				if(!is.matrix(X)) {
					X<-as.matrix(X)
				}
				if(nrow(X)!=S*R) {
					stop(paste("X must have S*R=",S*R," rows.",sep=""))
				}
				p<-ncol(X)
				if(qr(X)$rank!=p) {
					stop("Rank of X and its number of columns must be the same.")
				}
				if(qr(cbind(restr,X))$rank!=(S+p)) {
					stop("The columns of X are not linearly independent from the columns\nof the natural restrictions (kronecker(diag(S),rep(1,R))).")
				}
				XL<-A1%*%X
				if(!is.matrix(XL)) {
					XL<-as.matrix(XL)
				}
			} else {
				if(!is.matrix(XL)) {
					XL<-as.matrix(XL)
				}
				if(nrow(XL)!=u) {
					stop("XL and A1 must have the same number of rows.")
				}
				p<-ncol(XL)
				if(p>u) {
					stop("XL can't have more columns than rows.")
				}
				if(qr(XL)$rank!=p) {
					stop("Rank of XL and its number of columns must be the same.")
				}
			}
			UL<-t(COMPLEMENT(XL))
			form<-"freedom equation"

		} else {

			if(is.null(UL)) {
				if(!is.matrix(U)) {
					U<-t(U)
				}
				if(ncol(U)!=S*R) {
					stop("The number of columns of U have to be the same of the number of S*R.")
				}
				p<-S*(R-1)-nrow(U)
				if(p==0) {
					stop("The number of rows and columns of U must not be the same.")
				}
				if(qr(U)$rank!=(S*(R-1)-p)) {
					stop("Rank of U and its number of rows must be the same.")
				}
				XL<-A1%*%COMPLEMENT(rbind(U,t(restr)))
				UL<-t(COMPLEMENT(XL))
			} else {
				if(!is.matrix(UL)) {
					UL<-t(UL)
				}
				if(ncol(UL)!=u) {
					stop("The number of columns of UL have to be the same of the number of rows of A1.")
				}
				p<-u-nrow(UL)
				if(p==0) {
					stop("The number of rows and columns of UL must not be the same.")
				}
				if(qr(UL)$rank!=(u-p)) {
					stop("Rank of UL and its number of rows must be the same.")
				}
				XL<-COMPLEMENT(UL)
				if(!is.matrix(XL)) {
					XL<-as.matrix(XL)
				}
			}
			form<-"constraint"

		}
		if(all(XL==0)) {
			X<-0
		} else {
			X<-t(A1)%*%solve(A1%*%t(A1),tol=tole)%*%XL
		}
		Aa<-A1
		Xa<-X
		XLa<-XL
		if(u<S*(R-1)) {
			A0<-t(COMPLEMENT(rbind(A1,t(restr))))
			Aa<-rbind(A1,A0)
			if(all(XL==0)) {
				Xa<-cbind(t(A0)%*%solve(A0%*%t(A0),tol=tole))
				XLa<-rbind( matrix(0,u,S*(R-1)-u), diag(S*(R-1)-u) )
			} else {
				Xa<-cbind(X,t(A0)%*%solve(A0%*%t(A0),tol=tole))
				XLa<-rbind( cbind(XL,matrix(0,u,S*(R-1)-u)), cbind(matrix(0,S*(R-1)-u,p),diag(S*(R-1)-u)) )
			}
		}


	}

	if(type=="functional linear (obj)" | type=="functional linear (out)") {


		if(is.null(X) && is.null(U)) {
			stop("Specify X for freedom equation formulation or U for constraint formulation.")
		}
		if((!is.null(X)) && (!is.null(U))) {
			stop("Specify just X or U, not both.")
		}
		if(is.null(U)) {
			if(!is.matrix(X)) {
				X<-as.matrix(X)
			}
			if(nrow(X)!=u) {
				stop(paste("X must have the same number of rows of the length of the observed functions (=",u,").",sep=""))
			}
			p<-ncol(X)
			if(p>u) {
				stop("X can't have more columns than rows.")
			}
			if(qr(X)$rank!=p) {
				stop("Rank of X and its number of columns must be the same.")
			}
			U<-t(COMPLEMENT(X))
			form<-"freedom equation"
		} else {
			if(!is.matrix(U)) {
				U<-t(U)
			}
			if(ncol(U)!=u) {
				stop(paste("U must have the same number of columns of the length of the observed functions (=",u,").",sep=""))
			}
			p<-u-nrow(U)
			if(p==0) {
				stop("The number of rows and columns of U must not be the same.")
			}
			if(qr(U)$rank!=(u-p)) {
				stop("Rank of U and its number of rows must be the same.")
			}
			X<-COMPLEMENT(U)
			form<-"constraint"
		}


	}

	VFuI<-solve(VFu,tol=tole)
	glH<-u-p
	QwH<-0
	if(type=="linear") {
		if(u<S*(R-1)) {
			AA<-rbind(A,A0)
			XX<-rbind(cbind(X,matrix(0,u,S*(R-1)-u)),cbind(matrix(0,S*(R-1)-u,p),diag(S*(R-1)-u)))
		} else {
			AA<-A
			XX<-X
		}
		VthetaAI<-solve(AA%*%Vtheta%*%t(AA),tol=tole)
		Vbeta<-solve(t(XX)%*%VthetaAI%*%XX,tol=tole)
		beta<-c( Vbeta%*%t(XX)%*%VthetaAI%*%AA%*%theta )
		if(u<S*(R-1)){ thetabH<-Wpart%*%c(X%*%beta[1:p],rep(1,S),beta[(p+1):(S*(R-1)-u+p)]) } else { thetabH<-Wpart%*%c(X%*%beta,rep(1,S)) }
		thetaH<-obj$b+c( obj$B%*%thetabH )
		VthetabH<-W%*%Vbeta%*%t(W)
		VthetaH<-obj$B%*%VthetabH%*%t(obj$B)
		beta<-beta[1:p]
		Vbeta<-Vbeta[1:p,1:p]
		FH<-c(X%*%beta)  #F under GLM
		VFH<-X%*%Vbeta%*%t(X)
		if(glH>0) {
			QwH<-c(t(U%*%Fu)%*%solve(U%*%VFu%*%t(U),tol=tole)%*%U%*%Fu)
		}
	} else {
		if(type=="log-linear") {
			aux<-solve(Aa%*%solve(diag(theta))%*%Vtheta%*%solve(diag(theta))%*%t(Aa),tol=tole)
			Vbeta<-solve(t(XLa)%*%aux%*%XLa,tol=tole)
			beta<-c( Vbeta%*%t(XLa)%*%aux%*%Aa%*%log(theta) )
			thetaH<-c(solve(diag(c(kronecker(diag(S),matrix(1,R,R))%*%exp(Xa%*%beta))),tol=tole)%*%exp(Xa%*%beta))
			VLL<-matrix(0,S*R,S*R)
			for(s in 1:S) {
				VLL[((s-1)*R+1):(s*R),((s-1)*R+1):(s*R)]<-diag(thetaH[((s-1)*R+1):(s*R)])-thetaH[((s-1)*R+1):(s*R)]%*%t(thetaH[((s-1)*R+1):(s*R)])
			}
			VthetaH<-VLL%*%Xa%*%Vbeta%*%t(Xa)%*%VLL
			if(u<S*(R-1)) {
				if(p>0) {
					beta<-beta[1:p]
					Vbeta<-Vbeta[1:p,1:p]
				} else {
					beta<-0
					Vbeta<-0
				}
			}
			FH<-c(XL%*%beta)  #F under log-linear model
			VFH<-XL%*%Vbeta%*%t(XL)
			if(glH>0) {
				QwH<-c(t(UL%*%Fu)%*%solve(UL%*%VFu%*%t(UL),tol=tole)%*%UL%*%Fu)
			}
		} else {
			Vbeta<-solve(t(X)%*%VFuI%*%X,tol=tole)
			beta<-c( Vbeta%*%t(X)%*%VFuI%*%Fu )
			FH<-c(X%*%beta)  #F under the functional model
			VFH<-X%*%Vbeta%*%t(X)
			if(glH>0) {
				QwH<-c(t(U%*%Fu)%*%solve(U%*%VFu%*%t(U),tol=tole)%*%U%*%Fu)
			}
		}
	}

	if(type=="linear" | type=="log-linear") {
		ystH<-obj$Nst
		for(s in 1:S) {
			for(tt in 1:Tt[s]) {
				if(miss==TRUE && class(obj)=="satMarML" && obj$missing=="MAR") {
					if(tt==1) {
						zst<-diag(R)
					} else {
						zst<-obj$Zst[[paste("st",s,".",tt,sep="")]]
					}
					ystH[[paste("st",s,".",tt,sep="")]]<-c( obj$nsmm[s]*diag(thetaH[((s-1)*R+1):(s*R)])%*%zst%*%obj$alphast[[paste("st",s,".",tt,sep="")]] )
				} else {
					ystH[[paste("st",s,".",tt,sep="")]]<-obj$nstm[s,tt]*thetaH[((s-1)*R+1):(s*R)]
				}
			}
		}

		res<-list(obj=obj,R=R,S=S,miss=miss,Tt=Tt,zeroN=zeroN,theta=theta,Vtheta=Vtheta,call=match.call(),form=form,type=type,
			  p=p,u=u,thetaH=thetaH,VthetaH=VthetaH,beta=beta,Vbeta=Vbeta,Fu=Fu,VFu=VFu,FH=FH,VFH=VFH,QwH=QwH,glH=glH,ystH=ystH)
	} else {
		res<-list(obj=obj,zeroN=zeroN,theta=theta,Vtheta=Vtheta,call=match.call(),form=form,type=type,
			  p=p,u=u,beta=beta,Vbeta=Vbeta,Fu=Fu,VFu=VFu,FH=FH,VFH=VFH,QwH=QwH,glH=glH)
	}
	class(res)<-"funlinWLS"
	res
}


print.funlinWLS<-function(x, digits=max(3,getOption("digits")-3), ...) {
	cat("\nCall: ",deparse(x$call),"\n",sep="")

	if(x$form=="freedom equation") {
		pp<-matrix(0,x$p,4)
		pp[,1]<-x$beta
		pp[,2]<-sqrt(diag(as.matrix(x$Vbeta)))
		pp[,3]<-pp[,1]/pp[,2]
		pp[,4]<-1-pchisq((pp[,3])^2,1)
		colnames(pp)<-c("estimate","std.error","z-value","p-value")
		cat("\nWeighted least squares estimates of the parameters of the model:\n",sep="")
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)
	}

	cat("\nWald goodness of fit statistic of the model (d.f.=",x$glH,"): ",
	    round(x$QwH,digits=digits)," (p-value=",
	    ifelse(x$glH>0,format(round(1-pchisq(x$QwH,x$glH),digits=digits),digits=digits),1),")\n",sep="")

	cat("\n")
}


summary.funlinWLS<-function(object, digits=max(3,getOption("digits")-3), ...) {
	cat("\nCall: ",deparse(object$call),"\n",sep="")

	if(object$type=="linear" | object$type=="log-linear") {
		pp<-matrix(0,object$S,object$R)
		epp<-matrix(0,object$S,object$R)
		for(s in 1:object$S) {
			pp[s,]<-object$thetaH[((s-1)*object$R+1):(s*object$R)]
			epp[s,]<-sqrt(diag(object$VthetaH[((s-1)*object$R+1):(s*object$R),((s-1)*object$R+1):(s*object$R)]))
		}
		cat("\n\nWeighted least squares estimates of the probabilities under the model:\n",sep="")
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

		cat("\nStandard errors:\n",sep="")
		print.default(format(round(epp,digits=digits),digits=digits),quote=FALSE,print.gap=2)
	}

	pp<-matrix(0,object$u,4)
	pp[,1]<-object$Fu
	pp[,2]<-sqrt(diag(as.matrix(object$VFu)))
	pp[,3]<-object$FH
	pp[,4]<-sqrt(diag(as.matrix(object$VFH)))
	colnames(pp)<-c("observed","std.error","under the model","std.error")
	cat("\n\nWeighted least squares estimates of the functions:\n",sep="")
	print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)

	if(object$form=="freedom equation") {
		pp<-matrix(0,object$p,4)
		pp[,1]<-object$beta
		pp[,2]<-sqrt(diag(as.matrix(object$Vbeta)))
		pp[,3]<-pp[,1]/pp[,2]
		pp[,4]<-1-pchisq((pp[,3])^2,1)
		colnames(pp)<-c("estimate","std.error","z-value","p-value")
		cat("\n\nWeighted least squares estimates of the parameters of the model:\n",sep="")
		print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2)
	}

	cat("\n\nWald goodness of fit statistic of the model (d.f.=",object$glH,"): ",
	    round(object$QwH,digits=digits)," (p-value=",
	    ifelse(object$glH>0,format(round(1-pchisq(object$QwH,object$glH),digits=digits),digits=digits),1),")\n",sep="")

	if(object$type=="linear" | object$type=="log-linear") {
		if(object$miss==TRUE) {
			cat("\n\nAugmented estimated frequencies under the model:\n",sep="")
			for(s in 1:object$S) {
				if(object$S>1) {
					cat("\nSubpopulation ",s,"\n",sep="")
				}
				pp<-matrix(0,object$Tt[s],object$R)
				pp[1,]<-object$yst[[paste("st",s,".",1,sep="")]]
				if(object$Tt[s]>1) {
					for(tt in 2:object$Tt[s]) {
						pp[tt,]<-object$yst[[paste("st",s,".",tt,sep="")]]
					}
				}
				print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)
			}
		} else {
			cat("\n\nEstimated frequencies under the model:\n",sep="")
			pp<-matrix(0,object$S,object$R)
			for(s in 1:object$S) {
				pp[s,]<-object$yst[[paste("st",s,".",1,sep="")]]
			}
			print.default(format(round(pp,digits=digits),digits=digits),quote=FALSE,print.gap=2,right=TRUE)
		}
	}

	cat("\n")
}


#Wald statistic for testing a hypothesis H: C * beta = C0
waldTest<-function(obj,C,C0=NULL) {
	if(!any(class(obj)==c("linML","loglinML","funlinWLS"))) {
		stop("obj is not of the class linML, loglinML or funlinWLS")
	}
	if(obj$form!="freedom equation") {
		stop("obj must be a fitted model by freedom equation formulation")
	}
	tole<-1e-323 #tol argument for solve function

	p<-obj$p
	beta<-obj$beta
	Vbeta<-obj$Vbeta
	if(!is.matrix(C)) {
		C<-t(C)
	}
	if(ncol(C)!=p) {
		stop("The number of columns of C have to be the same of the number of parameters.")
	}
	nc<-nrow(C)
	if(nc>p) {
		stop("C can not have more rows than columns.")
	}
	if(qr(C)$rank!=nc) {
		stop("Rank of C and its number of rows must be the same.")
	}
	if(is.null(C0)) {
		C0<-numeric(nc)
	} else {
		if(is.vector(C0,mode="numeric")) {
			if(length(C0)!=nc) {
				stop("The length of C0 must be equal the number of rows of C")
			}
		} else {
			stop("C0 must be a numeric vector")
		}
	}
	QwH<-c(t(c(C%*%beta)-C0)%*%solve(C%*%Vbeta%*%t(C),tol=tole)%*%(c(C%*%beta)-C0))
	glH<-nc
	res<-list(call=match.call(),QwH=QwH,glH=glH)
	class(res)<-"waldTest"
	res
}

print.waldTest <- function(x,...){
				summary.waldTest(x,...)
}

summary.waldTest<-function(object, digits=max(3,getOption("digits")-3), ...) {
	cat("\nCall: ",deparse(object$call),"\n",sep="")

	cat("\nWald statistic of the hypothesis (d.f.=",object$glH,"): ",
	    round(object$QwH,digits=digits)," (p-value=",
	    ifelse(object$glH>0,format(round(1-pchisq(object$QwH,object$glH),digits=digits),digits=digits),1),")\n",sep="")

	cat("\n")
}


