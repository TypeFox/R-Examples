inQn<-function(x){
	fit2<-.C("R_inQn",as.integer(length(x)),as.double(x),as.double(0.0),DUP=TRUE,PACKAGE="DetMCD")
	fit2[[3]]
}
DetMCD<-function(X,h=NULL,alpha=0.75,scale_est="Auto",tol=1e-7){#h=NULL;alpha=0.5;scale_est="tau";tol=1e-7;X=x
	na.x<-complete.cases(X)
	if(sum(na.x)!=nrow(X)){	
		X<-X[na.x,]
		warning(paste('Observetions #:',which(na.x==0),"where removed"))	
	}
	Data<-data.matrix(X)
	n<-nrow(Data);
	p<-ncol(Data);
	if(n<(5*p))	stop("Too few observations")
	hf<-floor((n+p+1)/2)
	if(!is.numeric(alpha) & !is.numeric(h)) stop("alpha or h should be set")
	if(is.null(h)){
		alpha<-sort(alpha)
		if(min(alpha)>=0.5 & max(alpha)<=1){ 
			h<-sort(quanff(alpha,n=n,p=p))
		} else {
			stop("Error: invalid alpha value")
		}
	} 
	if(is.numeric(h)){
		h<-sort(h)
		if(min(h)<hf | max(h)>n) stop(paste("The smallest h should be at least ",hf," and at most ",n,sep=""))
	}
	if(min(h)<n){
		if(all(c("qn","tau","Auto")%in%scale_est==0)) stop("The scale_est not one of qn or tau or Auto")
		if(scale_est=="Auto") scale_est<-if(nrow(Data)>1000) "tau" else "qn"
		if(scale_est=="tau") scale_est<-"scaleTau2"
		if(scale_est=="qn") scale_est<-"qn"
		sca<-apply(Data,2,scale_est)	
		good_var<-(1:ncol(Data))[which(sca>tol)]
		if (length(good_var)<p) print(paste((1:p)[-good_var]," did not have positive variance",sep=""))
		if (length(good_var)<2) stop("'X' must have at least two columns with positive variance")
		Data<-Data[,good_var]	
		sca<-sca[good_var]
		med<-colMedians(Data)
		Data<-sweep(Data,2,med,FUN="-",check.margin=FALSE)
		Data<-sweep(Data,2,sca,FUN="/",check.margin=FALSE)
		out1<-DetMCD_SP(Data=Data,scale_est=scale_est,tol=tol)
		out2<-DetMCD_CS(Data=Data,scale_est=scale_est,h=h,out1=out1)
		out2$crit<-out2$objfun+2*sum(log(sca))
		out3<-lapply(1:length(h[h<n]),DetMCD_RW,Xw=X,scale_est=scale_est,out2=out2,alpha=alpha,hlst=h)
		if(max(h)==n) out3[[length(out3)+1]]<-DetMCD_RW(ll=1,hlst=min(h),Xw=X,out2=NULL,scale_est=scale_est,alpha=1)
	} else {
		out2<-vector("list",1)
		out2[[1]]<-DetMCD_RW(1,hlst=min(h),Xw=X,out2=NULL,scale_est=scale_est,alpha=1)
	}
	if(length(h)>1){
		outf<-xtractR_M(out3,X=X)
	} else {
		outf<-out3[[1]]
		outf$X<-X
	}
	class(outf)<-"DetMCD"
	outf
}
DetMCD_CS<-function(Data,scale_est,h,out1){
	n<-nrow(Data)
	p<-ncol(Data)
	hl<-length(h)
	hmin<-quanff(0.5,n,p);
	wM<-matrix(0,n,hl)
	if(scale_est=="qn") Mtype<-0 else Mtype<-1;
	fit2<-.C("R_FastR",
		as.integer(n),			#1
		as.integer(p),			#2
		as.double(Data),		#3
		as.double(out1),		#4
		as.integer(Mtype),		#5
		as.integer(h),			#6
		as.integer(hl),			#7
		as.integer(wM),			#8
		as.integer(rep(1,hl)),		#9
		as.integer(rep(1,hl)),		#10
		as.integer(ncol(out1)),		#11
		as.integer(rep(1,ncol(out1))),	#12
		as.double(rep(0.0,hl)),		#13
		as.integer(rep(0,hmin)),	#14
		as.integer(hmin),		#15
		as.integer(0),			#16
		DUP=TRUE,
	PACKAGE="DetMCD")
	if(fit2[[16]]==0){
		out<-list(mat=matrix(fit2[[8]],n,hl),which=fit2[[9]],nit=fit2[[10]],objfun=fit2[[13]],ef=fit2[[16]])
	} else {
		out<-list(mat=fit2[[14]],which=fit2[[9]][1],objfun=fit2[[13]],ef=fit2[[16]])
	}
	return(out)
}
quanff<-function(alpha,n,p) return(floor(2*floor((n+p+1)/2)-n+2*(n-floor((n+p+1)/2))*alpha))
DetMCD_SP<-function(Data,scale_est,tol){
	p<-ncol(Data)
	n<-nrow(Data)
	hsetsfull<-matrix(NA,p**2,6)
	#1: Htan of Data:
	hsetsfull[,1]<-c(cor(tanh(Data)))
	hsetsfull[,2]<-c(cor(Data,method="spearman"))
	#3: Tukey normal scores
	y3<-qnorm((apply(Data,2L,rank)-1/3)/(n+1/3))
	hsetsfull[,3]<-c(cor(y3,use="complete.obs"))
	#4: Spatial sign cov matrix
	znorm<-sqrt(rowSums(Data^2))
    	ii<-znorm>.Machine$double.eps
    	x.nrmd<-Data
    	x.nrmd[ii,]<-Data[ii,]/znorm[ii]
	hsetsfull[,4]<-c(crossprod(x.nrmd))
	#5: BACON
	ind5<-order(znorm);
    	half<-ceiling(n/2);
    	Hinit<-ind5[1:half]
	hsetsfull[,5]<-c(cov(Data[Hinit,,drop=FALSE]))
	Q<-rep(1.0,p**2)
	if(scale_est=="qn") Mtype<-0 else Mtype<-1;
	fit2<-.C("R_FastOGK",as.integer(n),as.integer(p),as.double(Data),as.double(Q),as.integer(Mtype),DUP=TRUE,PACKAGE="DetMCD")
	#matrix(fit2[[4]],p,p)
	#.25*(scaleTau2(Data[,1]+Data[,3])**2-scaleTau2(Data[,1]-Data[,3])**2);
	#l<-2;c<-5;.25*(qn(Data[,l]+Data[,c])**2-qn(Data[,l]-Data[,c])**2)/matrix(fit2[[4]],p,p)[l,c]
	hsetsfull[,6]<-fit2[[4]]
	return(hsetsfull)
}
DetMCD_tst<-function(Data,scale_est,tol){
	p<-ncol(Data)
	n<-nrow(Data)
	Q<-rep(1.0,p**2)
	if(scale_est=="qn") Mtype<-0 else Mtype<-1;
	fit2<-.C("R_FastTEST",as.integer(n),as.integer(p),as.double(Data),as.double(Q),as.integer(Mtype),PACKAGE="DetMCD")
	#matrix(fit2[[4]],p,p)
	#.25*(scaleTau2(Data[,1]+Data[,2])**2-scaleTau2(Data[,1]-Data[,2])**2);
	#.25*(qn(Data[,1]+Data[,2])**2-qn(Data[,1]-Data[,2])**2)
	hsetsfull[,6]<-fit2[[4]]
	return(hsetsfull)
}
DetMCD_RW<-function(ll,hlst,Xw,out2=NULL,scale_est,alpha){	
	p<-ncol(Xw);
	n<-nrow(Xw)
	if(is.null(out2))	out2$ef<-0
	if(out2$ef==0){
		h<-hlst[ll]
		if(!is.null(out2$mat))	Isets<-out2$mat[1:h,ll] else Isets<-1:n
		mah<-mahalanobis(Xw,colMeans(Xw[Isets,]),var(Xw[Isets,]))	
		factor<-quantile(mah,min(h/n,(h-1)/n),type=1)/qchisq(min(h/n,(h-1)/n),p)
		raw.cov<-factor*var(Xw[Isets,])
		raw.center<-colMeans(Xw[Isets,])
		raw.objective<-log(det(var(Xw[Isets,])))
		mah<-mah/factor
		raw.rd<-sqrt(mah)
		cutoff.rd<-sqrt(qchisq(0.975,df=p))
		weights<-as.numeric(raw.rd<=cutoff.rd)
		raw.wt<-weights
		wstats<-cov.wt(Xw,wt=weights,center=TRUE)
		rew.center<-wstats$center
		rew.cov<-wstats$cov
		mah<-mahalanobis(Xw,rew.center,rew.cov)
		rew.rd<-sqrt(mah)
		rew.flag<-as.numeric(rew.rd<=cutoff.rd)
		names(raw.center)<-names(Xw)		
		dimnames(raw.cov)<-list(dimnames(Xw)[[2]],dimnames(Xw)[[2]])
		HrewS<-(1:n)[which(raw.wt>0)]
		if(!is.null(out2$mat))	best_sub<-out2$which else best_sub<-1
		out<-list(center=rew.center,cov=rew.cov,Hsubsets=HrewS,rd=rew.rd,best_sub=best_sub,
		flag=rew.flag,raw.center=raw.center,raw.cov=raw.cov,best=Isets,cutoff=cutoff.rd,
		crit=raw.objective,h=h,raw.rd=raw.rd,raw.wt=raw.wt,scale_est=scale_est,alpha=alpha[ll])
	} else {
		warning("More than [(n+p+1)/2] of the observations lie on a hyperplane.")
		Isets<-out2$mat
		out<-list(center=colMeans(Xw[Isets,]),crit=-Inf,cov=var(Xw[Isets,]),Hsubsets=Isets,best_sub=out2$which)
	}
	return(out)
}
xtractR_M<-function(out2,X){
	nEntry<-length(out2)
	nH<-length(out2[[1]]$Hsubsets)
	p<-length(out2[[1]]$center)
	cutoff<-slnames<-rep(NA,nEntry)
	trcenter<-raw.center<-matrix(NA,nEntry,p)
	rew.rd<-rew.flag<-raw.rd<-raw.wt<-matrix(NA,nrow(X),nEntry)
	best<-matrix(NA,nrow(X),nEntry)
	rn<-row.names(out2[[1]]$raw.cov)
	raw.cov<-trcov<-array(NA,c(dim(out2[[1]]$cov),nEntry))
	alpha<-raw.objective<-scale_est<-rep(NA,nEntry)
	hlist<-rep(NA,nEntry);
	for(i in 1:nEntry){
		alpha[i]<-out2[[i]]$alpha
		slnames[i]<-paste("h_",i,sep="")
		trcenter[i,]<-out2[[i]]$center
		trcov[,,i]<-out2[[i]]$cov
		lh<-length(out2[[i]]$best)
		best[1:lh,i]<-out2[[i]]$best
		rew.rd[,i]<-out2[[i]]$rd
		rew.flag[,i]<-out2[[i]]$flag
		raw.center[i,]<-out2[[i]]$raw.center
		raw.cov[,,i]<-out2[[i]]$raw.cov
		raw.objective[i]<-out2[[i]]$crit
		hlist[i]<-out2[[i]]$h
		scale_est[i]<-out2[[i]]$scale_est
		raw.rd[,i]<-out2[[i]]$raw.rd
		raw.wt[,i]<-out2[[i]]$raw.wt
		cutoff[i]<-out2[[i]]$cutoff
	}
	colnames(trcenter)<-colnames(raw.center)<-names(out2[[i]]$center)
	rownames(trcenter)<-rownames(raw.center)<-slnames
	colnames(best)<-slnames
	colnames(rew.flag)<-colnames(rew.rd)<-slnames
	colnames(raw.wt)<-colnames(raw.rd)<-slnames
	names(raw.objective)<-names(hlist)<-names(scale_est)<-slnames
	if(!is.null(alpha))	names(alpha)<-slnames
	dimnames(trcov)<-dimnames(raw.cov)<-list(rn,rn,slnames)
	out3<-list(raw.center=raw.center,raw.cov=raw.cov,crit=raw.objective,
	raw.rd=raw.rd,raw.wt=t(raw.wt),center=trcenter,cov=trcov,h=hlist,rd=rew.rd,
	weights=t(rew.flag),scale_est=scale_est,X=X,alpha=alpha,best=best)
	return(out3)
}
plot.DetMCD<-function(x,h.val=1,which=c("all","dd","distance","qqchi2","tolEllipsePlot","screeplot"),classic=FALSE,
ask=(which=="all"&&dev.interactive()),cutoff=NULL,id.n,labels.id=rownames(x),cex.id=0.75,label.pos=c(4,2),tol=1e-7,...){
#Original code from robustbase::covPlot. See citation("robustbase");
	h.val<-as.integer(h.val)
	h.max<-nrow(x[[1]])
	if(!is.null(h.max)){
		if(is.numeric(h.val)==0 | h.val>h.max){ 
			stop(paste("h.val must be a integer smaller or equal to ",h.max,sep=""))
		}
		m.cov=list(center=x$center[h.val,],cov=x$cov[,,h.val])
	} else {
		m.cov=list(center=x$center,cov=x$cov)
	}
	Data.X<-x$X
	myscreeplot<-function(Data.X,m.cov){
		erob<-eigen(m.cov$cov,symmetric=TRUE,only.values=TRUE)$values
		eclass<-eigen(var(Data.X),symmetric=TRUE,only.values=TRUE)$values
		leg.txt<-c("Robust","Classical")
		leg.col<-c("green","red")
		leg.pch<-c(1,24)
		leg.lty<-c("solid","dotted")
		eall<-c(erob,eclass)
		ylim<-c(min(eall),max(eall))
		plot(erob,ylim=ylim,ylab="Eigenvalues",xlab="Index",type="n")
		legend("topright",leg.txt,pch=leg.pch,lty=leg.lty,col=leg.col)
		lines(erob,type="o",pch=leg.pch[1],lty=leg.lty[1],col=leg.col[1])
		lines(eclass,type="o",pch=leg.pch[2],lty=leg.lty[2],col=leg.col[2])
		title(main="Scree plot")
	}
	mydistplot<-function(x,cutoff,classic=FALSE,id.n){
		n<-length(x)
		if(missing(id.n)) 
		id.n<-length(which(x > cutoff))
		ylab<-paste("Square Root of",if(classic) "Mahalanobis" else "Robust","distance")
		plot(x,type="p",ylab=ylab,xlab="Index",main="Distance Plot")
		label(1:n,x,id.n)
		abline(h=cutoff)
    	}
	myddplot<-function(md,rd,cutoff,id.n){
		n<-length(md)
		if(missing(id.n)) id.n<-length(which(rd>cutoff))
		xlab<-"Mahalanobis distance"
		ylab<-"Robust distance"
		plot(md,rd,type="p",xlab=xlab,ylab=ylab,main="Distance-Distance Plot")
		label(md,rd,id.n)
		abline(0,1,lty=2)
		abline(v=cutoff,h=cutoff)
    	}
	qqplot<-function(x,p,cutoff=sqrt(qchisq(0.975,p)),classic=FALSE,id.n){
		n<-length(x)
		if(missing(id.n)) id.n<-length(which(x > cutoff))
		qq<-sqrt(qchisq(((1:n)-1/3)/(n+1/3),p))
		x<-sort(x,index.return=TRUE)
		ix<-x$ix
		x<-x$x
		ylab<-paste(if(classic) "Mahalanobis" else "Robust","distance")
		xlab<-"Square root of the quantiles of the chi-squared distribution"
		plot(qq,x,xlab=xlab,ylab=ylab,main="Chisquare QQ-Plot")
		label(qq,x,id.n,ind=(n - id.n + 1):n,labs=ix)
		abline(0,1,lty=2)
	}
	label<-function(x,y,id.n,ind=sort.list(y,decreasing=TRUE)[1:id.n],labs=labels.id,adj.x=TRUE){
		if(id.n>0){
			labpos<-if(adj.x) label.pos[1 + as.numeric(x > mean(range(x)))]	else 3
		    	text(x[ind],y[ind],labs[ind],cex=cex.id,xpd=TRUE,pos=labpos,offset=0.25)
		}
    	}
	if(is.data.frame(Data.X)) Data.X<-data.matrix(Data.X)
	if(!is.matrix(Data.X) || !is.numeric(Data.X)) stop("x is not a numeric Dataframe or matrix.")
	n<-dim(Data.X)[1]
	p<-dim(Data.X)[2]
	if(!is.numeric(m.cov$center) || !is.numeric(m.cov$cov)){ 
		stop("argument 'm.cov' must have numeric components 'center' and 'cov'")
	}
	if(length(m.cov$center)!=p) stop("Data set and provided center have different dimensions!")
    	if(is.numeric(m.cov$crit) && m.cov$crit==0) stop("The covariance matrix is singular!")
    	if(is.null(cutoff)) cutoff<-sqrt(qchisq(0.975,p))
    	if(is.null(labels.id)) labels.id<-as.character(1:n)
    	if(!missing(id.n) && !is.null(id.n)){
		id.n<-as.integer(id.n)
		if(id.n<0 || id.n>n) stop(sQuote("id.n")," must be in{1,..,",n,"}")
    	}
    	which<-match.arg(which)
    	md<-sqrt(mahalanobis(Data.X,colMeans(Data.X),var(Data.X)))
    	rd<-sqrt(mahalanobis(Data.X,m.cov$center,m.cov$cov))
    	op<-if(ask) par(ask=TRUE) else list()
    	on.exit(par(op))
    	if(which == "all" || which == "distance"){
		if(classic){
    			opr<-if(prod(par("mfrow"))==1) 
			par(mfrow=c(1,2),pty="m")
    		else list()
		}
		mydistplot(rd,cutoff,id.n=id.n)
		if(classic){
		    	mydistplot(md,cutoff,classic=TRUE,id.n=id.n)
		    	par(opr)
		}
    	}
    	if(which=="all" || which=="dd") myddplot(md,rd,cutoff=cutoff,id.n=id.n)
    	if(which=="all" || which=="qqchi2"){
	if(classic){
    		opr<-if(prod(par("mfrow"))==1) 
			par(mfrow=c(1,2),pty="m")
    		else 	list()
	}
	qqplot(rd,p,cutoff=cutoff,id.n=id.n)
	if(classic){
    		qqplot(md,p,cutoff=cutoff,classic=TRUE,id.n=id.n)
    		par(opr)
	}
   	}
    	if(which=="all" || which=="tolEllipsePlot"){
		if(p==2){
			robustbase::tolEllipsePlot(Data.X,m.cov=m.cov,cutoff=cutoff,id.n=id.n,classic=classic,tol=tol)
		} else 
			if(which!="all") warning("For tolerance ellipses the dimension 'p' must be 2!")
    	}
    	if(which=="all" || which=="screeplot"){
		myscreeplot(Data.X,m.cov=m.cov)
    	}
}
