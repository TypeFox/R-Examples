library(impute)
library(Biobase)
library(combinat)
MetaDE.rawdata<-function(x, ind.method=c("modt","regt","pairedt","F","pearsonr","spearmanr","logrank"),
                    meta.method=c("maxP","maxP.OC","minP","minP.OC","Fisher","Fisher.OC","AW","AW.OC"
                                  ,"roP","roP.OC","Stouffer","Stouffer.OC","SR","PR","minMCC","FEM","REM","rankProd"),
                    paired=NULL,miss.tol=0.3,rth=NULL,nperm=NULL,ind.tail="abs",asymptotic=FALSE,...)
{
	x<-MetaDE.impute(x,y=miss.tol)
   ind.method<-match.arg(ind.method,c("regt","modt","pairedt","pearsonr","spearmanr","F","logrank"),several.ok=TRUE)
   meta.method<-match.arg(meta.method,c("maxP","maxP.OC","minP","minP.OC","Fisher","Fisher.OC","AW","AW.OC"
                                  ,"roP","roP.OC","Stouffer","Stouffer.OC","SR","PR","minMCC","FEM","REM","rankProd"),several.ok = TRUE)
   check.dim(x,ind.method=ind.method,meta.method=meta.method,paired=paired)
   check.method1(x,ind.method=ind.method,meta.method=meta.method,rth=rth,paired=paired)
   check.asymptotic(meta.method=meta.method,asymptotic=asymptotic)
   check.tail(meta.method=meta.method,ind.tail)
   x<-check.exp(x) 
   if ("minMCC"%in%meta.method){
		K<-length(x)
		dat<-lbl<-list()
		for(i in 1:K){
			dat[[i]]<-x[[i]][[1]]
			lbl[[i]]<-x[[i]][[2]]
		}
		res<-get.minMCC(dat=data,lbl=lbl,nperm=nperm)
    }else if ("rankProd"%in%meta.method){
		K<-length(x)
		dat<-lbl<-list()
		for(i in 1:K){
			dat[[i]]<-x[[i]][[1]]
			lbl[[i]]<-x[[i]][[2]]
		}
		res<-get.RP(dat=dat,lbl=lbl,nperm=nperm)
   }else if ("FEM"%in%meta.method|"REM"%in%meta.method){
	ind.res<-ind.cal.ES(x,paired=paired,nperm=nperm,miss.tol=miss.tol)
	meta.res<-MetaDE.ES(ind.res,meta.method=meta.method)
	res<-list(meta.analysis=meta.res,ind.ES=ind.res$ES,ind.Var=ind.res$Var,ind.stat=NULL,ind.p=NULL,raw.data=x)
	attr(res$meta.analysis,"nperstudy")<-attr(ind.res,"nperstudy")
    attr(res$meta.analysis,"nperlabelperstudy")<-attr(ind.res,"nperlabelperstudy")
   }else {
   ind.res<-ind.analysis(x,ind.method=ind.method,nperm=nperm,miss.tol=miss.tol,tail=ind.tail,...)   
  if (asymptotic) 
   {
    ind.res$bp<-NULL
    cat("Asymptotic estimation was used instead of the permutation\n")
   }else cat("Permutation was used instead of the asymptotic estimation\n")
    nm<-length(meta.method)
    meta.res<-MetaDE.pvalue(ind.res,meta.method=meta.method,asymptotic=asymptotic,rth=rth,miss.tol=miss.tol)$meta.analysis
    res<-list(meta.analysis=meta.res,ind.stat=ind.res$stat,ind.p=ind.res$p,raw.data=x)
  }
  if("minMCC"%in%meta.method){
		class(res)<-"MetaDE.minMCC"
	}else if("FEM"%in%meta.method|"REM"%in%meta.method){
		class(res)<-"MetaDE.ES"
	}else{
		class(res)<-"MetaDE.pvalue"
	}
   return(res)
}
  
MetaDE.minMCC<-function(x,nperm=100,miss.tol=0.3)
{
   x<-check.exp(x)
   x<-MetaDE.impute(x,y=miss.tol)
   K<-length(x)
   dat<-lbl<-list()
   N<-n<-NULL
   for(i in 1:K){
		dat[[i]]<-x[[i]][[1]]
		lbl[[i]]<-x[[i]][[2]]
            nns<-get.sample.label.number2(lbl[[i]])
            N<-c(N,nns$N)
            n<-c(n,nns$n)
	}
   meta.res<-get.minMCC(dat,lbl,nperm=nperm)
   colnames(meta.res$stat)<-colnames(meta.res$pval)<-colnames(meta.res$FDR)<-"minMCC"
   attr(meta.res,"nperstudy")<-N
   attr(meta.res,"nperlabelperstudy")<-n 
   res<-list(meta.analysis=meta.res,raw.data=x)
   attr(res$meta.analysis,"nperstudy")<-N
   res$ind.stat<-NA
   res$ind.p<-NA
   class(res)<-"MetaDE.minMCC"
   return(res)
}
   
MetaDE.pvalue <-function(x,meta.method=c("maxP","maxP.OC","minP","minP.OC","Fisher","Fisher.OC","AW","AW.OC"
         ,"roP","roP.OC","Stouffer","Stouffer.OC","SR","PR"),rth=NULL,miss.tol=0.3,asymptotic=FALSE)
{
     meta.method<-match.arg(meta.method,several.ok = TRUE)
     check.asymptotic(meta.method,asymptotic)
      K<-ncol(x$p)
     if (asymptotic) x$bp<-NULL
       
      nm<-length(meta.method)
      meta.res<-list(stat=NA,pval=NA,FDR=NA,AW.weight=NA)
      meta.res$stat<-meta.res$pval<-meta.res$FDR<-matrix(NA,nrow(x$p),nm)
     for( i in 1:nm){
       temp<-switch(meta.method[i],maxP={get.maxP(x$p,x$bp,miss.tol=miss.tol)},minP={get.minP(x$p,x$bp,miss.tol=miss.tol)},
                      Fisher={get.fisher(x$p,x$bp,miss.tol=miss.tol)},roP={get.roP(x$p,x$bp,rth=rth)},
                      AW={get.AW(x$p,x$bp)},AW.OC={get.AW.OC(x$p,x$bp)},Fisher.OC={get.fisher.OC(x$p,x$bp)},
                      maxP.OC={get.maxP.OC(x$p,x$bp,miss.tol=miss.tol)},minP.OC={get.minP.OC(x$p,x$bp,miss.tol=miss.tol)},
                      roP.OC={get.roP.OC(x$p,x$bp,rth=rth)},
                      Stouffer={get.Stouff(x$p,x$bp,miss.tol=miss.tol)},Stouffer.OC={get.Stouff.OC(x$p,x$bp)},
					  SR={get.SR(x$p,x$bp)},PR={get.PR(x$p,x$bp)})
       meta.res$stat[,i]<-temp$stat
       meta.res$pval[,i]<-temp$pval
       meta.res$FDR[,i]<-temp$FDR
	   if(meta.method[i]=="AW"|meta.method[i]=="AW.OC"){meta.res$AW.weight<-temp$AW.weight}
      }
     colnames(meta.res$stat)<-colnames(meta.res$pval)<-colnames(meta.res$FDR)<-meta.method
     rownames(meta.res$stat)<-rownames(meta.res$pval)<-rownames(meta.res$FDR)<-rownames(x$p)   
	 attr(meta.res,"nstudy")<-K
     attr(meta.res,"meta.method")<-meta.method 
	res<-list(meta.analysis=meta.res,ind.p=x$p)	 
	class(res)<-"MetaDE.pvalue"
     return(res)
}
 MetaDE.ES<-function(x,meta.method=c("FEM","REM")){
	meta.method<-match.arg(meta.method)
	K<-ncol(x$ES)
	if(meta.method=="REM"){
		res<-get.REM(x$ES,x$Var,pe=x$perm.ES,pv=x$perm.Var)
		tempFDR<-matrix(res$FDR,ncol=1)
		rownames(tempFDR)<-rownames(x$ES)
		colnames(tempFDR)<-meta.method
		meta.res<-list(mu.hat=res$mu.hat,mu.var=res$mu.var,Qval=res$Qval,Qpval=res$Qpval,tau2=res$tau2,zval=res$zval,pval=res$pval,FDR=tempFDR)	
	}else{
		res<-get.FEM(x$ES,x$Var,x$perm.ES,x$perm.Var)
		tempFDR<-matrix(res$FDR,ncol=1)
		rownames(tempFDR)<-rownames(x$ES)
		colnames(tempFDR)<-meta.method
		meta.res<-list(mu.hat=res$mu.hat,mu.var=res$mu.var,zval=res$zval,pval=res$pval,FDR=tempFDR)	
	}	
    attr(meta.res,"nstudy")<-K
    attr(meta.res,"meta.method")<-meta.method 
	class(meta.res)<-"MetaDE.ES"
    return(meta.res)
}
get.fisher<-function(p,bp=NULL,miss.tol) {
    k<-ncol(p)
	pval<-stat<-rep(NA,nrow(p))
	if(!is.null(bp)){
		rnum<-which(apply(p,1,function(x) !any(is.na(x))))
		Ug<-matrix(NA,nrow(p),1)
		Ug[rnum,1]<-as.matrix(rowSums(-2*log(p[rnum,])))
		Ubg<-matrix(rowSums(-2*log(bp)),nrow(p),nrow(bp)/nrow(p))
		pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="high")
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=Ug, pval=pval, FDR=qval)
	}else{
		rnum<-which(apply(p,1,function(x) sum(is.na(x))/k)<miss.tol)
		pval[rnum]<-apply(p[rnum,],1,function(x)1-pchisq(-2*sum(log(x),na.rm=T),2*sum(!is.na(x))))
		qval<-p.adjust(pval,method="BH")
		stat[rnum]<-apply(p[rnum,],1,function(x)-2*sum(log(x),na.rm=T))
		res<-list(stat=stat,pval=pval,FDR=qval)
	}
	names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
    return(res)
}
get.fisher.OC<-function(p,bp=NULL)
{
     K<-ncol(p)
	 rnum<-which(apply(p,1,function(x) !any(is.na(x))))
	 pval<-rep(NA,nrow(p))
	 Ug<-matrix(NA,nrow(p),1)
	 Ug[rnum,1]<-pmax(rowSums(-log(p))[rnum],rowSums(-log(1-p))[rnum])
      if (is.null(bp)) 
      {
       warning("there're no parametric results for Pearson's method,we will use simulation to estimate 
       the p values")
       bp<-matrix(runif(500*nrow(p)*K),500*nrow(p),K) 
      }
        Ubg<-matrix(pmax(rowSums(-log(bp)),rowSums(-log(1-bp))),nrow(p),nrow(bp)/nrow(p))
	 pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="high")
        qval<-p.adjust(pval,method="BH")

      res<-list(stat=Ug,pval=pval,FDR=qval)
	 names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
      return(res)
}

get.maxP<-function(p,bp=NULL,miss.tol) {
	k<-ncol(p)
	pval<-stat<-rep(NA,nrow(p))
	if(!is.null(bp)){
		rnum<-which(apply(p,1,function(x) !any(is.na(x))))
		Ug<-matrix(NA,nrow(p),1)
		Ug[rnum,1]<-as.matrix(apply(p[rnum,],1,max))
		Ubg<-matrix(apply(bp,1,max),nrow(p),nrow(bp)/nrow(p))
		pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="low")
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=Ug, pval=pval, FDR=qval)
	}else{
		rnum<-which(apply(p,1,function(x) sum(is.na(x))/k)<=miss.tol)
		pval[rnum]<-apply(p[rnum,],1,function(x)pbeta(max(x,na.rm=T),sum(!is.na(x)),1))
		stat[rnum]<-apply(p[rnum,],1,function(x)max(x,na.rm=T))
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=stat,pval=pval,FDR=qval)
	}
	names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
    return(res)
}


get.maxP.OC<-function(p,bp=NULL,miss.tol) {
	k<-ncol(p)
	pval<-stat<-rep(NA,nrow(p))
	if(!is.null(bp)){
		rnum<-which(apply(p,1,function(x) !any(is.na(x))))
		Ug<-matrix(NA,nrow(p),1)
		Ug[rnum,1]<-apply(p[rnum,],1,function(x)min(sort(x)[k],sort(1-x)[k]))
		Ubg<-matrix(apply(bp,1,function(x)min(sort(x)[k],sort(1-x)[k])),nrow(p),nrow(bp)/nrow(p))
		pval[rnum]<-perm.p(Ug[rnum,],Ubg[rnum,],tail="low")
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=c(Ug), pval=pval, FDR=qval)
	}else{
		rnum<-which(apply(p,1,function(x) sum(is.na(x))/k)<=miss.tol)
		stat[rnum]<-apply(p[rnum,],1,function(x)min(sort(x)[sum(!is.na(x))],sort(1-x)[sum(!is.na(1-x))]))
		k0<-apply(p[rnum,],1,function(x)sum(!is.na(x)))
		pval[rnum]<-ifelse(stat[rnum]>0.5,2*stat[rnum]^k0-(2*stat[rnum]-1)^k0,2*stat[rnum]^k0)	
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=stat,pval=pval,FDR=qval)
	}
	names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
    return(res)
}


get.minP<-function(p,bp=NULL,miss.tol){
    k<-ncol(p)
	pval<-stat<-rep(NA,nrow(p))
    if (!is.null(bp)){
		rnum<-which(apply(p,1,function(x) !any(is.na(x))))
		Ug<-matrix(NA,nrow(p),1)
		Ug[rnum,1]<-apply(p[rnum,],1,min)
		Ubg<-matrix(apply(bp,1,min),nrow(p),nrow(bp)/nrow(p))
        pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="low")
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=c(Ug), pval=pval, FDR=qval)
      }else{
		rnum<-which(apply(p,1,function(x) sum(is.na(x))/k)<=miss.tol)
		pval[rnum]<-apply(p[rnum,],1,function(x)pbeta(min(x,na.rm=T),1,sum(!is.na(x))))
		stat[rnum]<-apply(p[rnum,],1,function(x)min(x,na.rm=T))
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=stat,pval=pval,FDR=qval)
	}
	names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
    return(res)
}


get.minP.OC<-function(p,bp=NULL,miss.tol)
{
    K<-ncol(p)
	pval<-stat<-rep(NA,nrow(p))
    if (!is.null(bp)){
		rnum<-which(apply(p,1,function(x) !any(is.na(x))))
		Ug<-matrix(NA,nrow(p),1)
		Ug[rnum,1]<-as.matrix(apply(cbind(p,1-p)[rnum,],1,min))
		Ubg<-matrix(apply(cbind(bp,1-bp),1,min),nrow(p),nrow(bp)/nrow(p))
        pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="low")
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=c(Ug),pval=pval,FDR=qval)
      }else{
	  	rnum<-which(apply(p,1,function(x) sum(is.na(x))/K)<=miss.tol)
	  	stat[rnum]<-apply(cbind(p,1-p)[rnum,],1,min,na.rm=T)
		K0<-apply(p[rnum,],1,function(x)sum(!is.na(x)))
		pval[rnum]<-ifelse(stat[rnum]>=0.5,1,1-(1-2*stat[rnum])^K0)
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=stat,pval=pval,FDR=qval)
      }
      names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
      return(res)
}

get.roP<-function(p,bp=NULL,rth) {
	k<-ncol(p)
	rnum<-which(apply(p,1,function(x) !any(is.na(x))))
	pval<-stat<-rep(NA,nrow(p))
	if(!is.null(bp)){
		p<-t(apply(p,1, sort,na.last = T))
		bp<-t(apply(bp,1,sort,na.last = T))
		Ug<-matrix(NA,nrow(p),1)
		Ug[rnum,1]<-p[rnum,rth]	
		Ubg<-matrix(bp[,rth],nrow(p),nrow(bp)/nrow(p))
		pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="low")
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=Ug,pval=pval, FDR=qval)
	}else{
		pval[rnum]<-apply(p[rnum,],1,function(x)pbeta(x[order(x)][rth],rth,k-rth+1))			
		stat[rnum]<-apply(p[rnum,],1,function(x)x[order(x)][rth])			
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=stat,pval=pval,FDR=qval)	
	}
	names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
    return(res)
}

get.roP.OC<-function(p,bp=NULL,rth) {
	k<-ncol(p)
	rnum<-which(apply(p,1,function(x) !any(is.na(x))))
	pval<-stat<-rep(NA,nrow(p))
	if(!is.null(bp)){
		Ug<-matrix(NA,nrow(p),1)
		Ug[rnum,1]<-apply(p[rnum,],1,function(x)min(sort(x)[rth],sort(1-x)[rth]))
		Ubg<-matrix(apply(bp,1,function(x)min(sort(x)[rth],sort(1-x)[rth])),nrow(p),nrow(bp)/nrow(p))
		pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="low")
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=c(Ug), pval=pval, FDR=qval)
	}else{
		stat[rnum]<-apply(p[rnum,],1,function(x)min(sort(x)[rth],sort(1-x)[rth]))
		pval[rnum]<-sapply(stat[rnum],function(x)CDF.rop(x,rth,k))	
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=stat,pval=pval,FDR=qval)
	}
	names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
    return(res)
}


get.Stouff<-function(p,bp=NULL,miss.tol){
	k<-ncol(p)
	pval<-stat<-rep(NA,ncol(p))
	if(!is.null(bp)){
		rnum<-which(apply(p,1,function(x) !any(is.na(x))))
		Ug<-matrix(NA,nrow(p),1)
		Ug[rnum,1]<-apply(p[rnum,],1,function(x)sum(qnorm(x))/sqrt(k))
		Ubg<-matrix(apply(bp,1,function(x)sum(qnorm(x))/sqrt(k)),nrow(p),nrow(bp)/nrow(p))
		pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="abs")
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=c(Ug), pval=pval, FDR=qval)
	}else{
		rnum<-which(apply(p,1,function(x) sum(is.na(x))/k)<=miss.tol)
		pval[rnum]<-apply(p[rnum,],1,function(x)2*(1-pnorm(abs(sum(qnorm(x),na.rm=T)/sqrt(sum(!is.na(x)))))))
		stat[rnum]<-apply(p[rnum,],1,function(x)sum(qnorm(x),na.rm=T)/sqrt(sum(!is.na(x))))
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=stat,pval=pval,FDR=qval)
	}
	names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
    return(res)
}


get.Stouff.OC<-function(p,bp=NULL){
	k<-ncol(p)
	rnum<-which(apply(p,1,function(x) !any(is.na(x))))
    pval<-rep(NA,nrow(p),1)
	Ug<-UL<-UR<-matrix(NA,nrow(p),1)
	UL[rnum,1]<-apply(p[rnum,],1,function(x)sum(qnorm(x))/sqrt(k))
	UR[rnum,1]<-apply((1-p)[rnum,],1,function(x)sum(qnorm(x))/sqrt(k))
	Ug[rnum,1]<-pmax(UL,UR)[rnum]
	if(!is.null(bp)){
		UbL<-as.matrix(apply(bp,1,function(x)sum(qnorm(x))/sqrt(k)))
		UbR<-as.matrix(apply(1-bp,1,function(x)sum(qnorm(x))/sqrt(k)))
		Ubg<-pmax(UbL,UbR)
		pval[rnum]<-perm.p(Ug[rnum,1],Ubg[rnum,],tail="high")
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=c(Ug), pval=pval, FDR=qval)
	}else{
		pval[rnum]<-sapply(Ug[rnum,1],function(x)ifelse(x>0,2*(1-pnorm(x)),0))
		qval<-p.adjust(pval,method="BH")
		res<-list(stat=c(Ug),pval=pval,FDR=qval)
	}
	names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
    return(res)
}

get.AW<-function(p,bp=NULL) {
    gene.names<-row.names(p) # get gene names
    K<-ncol(p)
    rnum<-which(apply(p,1,function(x) !any(is.na(x))))

      K<-ncol(p)
      if (is.null(bp)) 
      {
       warning("there're no parametric results for Pearson's method,we will use simulation to estimate 
       the p values")
       bp<-matrix(runif(5000*K),5000,K) 
      }

    nW<-2^K-1
    W<-do.call(expand.grid, rep(list(c(0, 1)), K))[-1,]
    V<--log(p)%*%t(W)
    Vb<--log(bp)%*%t(W)
    nr<-nrow(V)
    np<-nrow(Vb)
    minV<-index.AW<-pval<-qval<-matrix(NA,nr,1)
    pV<-matrix(NA,nr,nW)
    minVb<-matrix(9999,nr,np/nr) #give a large initial value
    minVb[-rnum,]<-NA
    for (i in 1:nW) {
        m.Vb<-matrix(Vb[,i],nr,np/nr)
        pV[rnum,i]<-perm.p(as.matrix(V[rnum,i]),m.Vb[rnum,],tail="high")
        minVb<-pmin(minVb ,emp(m.Vb[rnum,],tail="high"))
    }
    minV[rnum]<-apply(pV[rnum,],1,min,na.rm=T)
    index.AW[rnum]<-apply(pV[rnum,],1,which.min)
    AW.weight<-W[index.AW,]
    pval[rnum]<-perm.p(minV[rnum],minVb[rnum,],tail="low")
    qval[rnum]<-p.adjust(pval[rnum],method="BH")
    row.names(minV)<-row.names(pval)<-row.names(qval)<-row.names(AW.weight)<-gene.names
    res<-list(stat=minV, pval=pval, FDR=qval, AW.weight=AW.weight)
    return(res)
}



get.AW.OC<-function(p,bp=NULL) {
    gene.names<-row.names(p) 
    K<-ncol(p)
    rnum<-which(apply(p,1,function(x) !any(is.na(x))))
      K<-ncol(p)
      if (is.null(bp)) 
      {
       warning("there're no parametric results for Pearson's method,we will use simulation to estimate 
       the p values")
       bp<-matrix(runif(5000*K),5000,K) 
      }

    nW<-2^K-1
	W<-do.call(expand.grid, rep(list(c(0, 1)), K))[-1,]
	Vl<--log(p)%*%t(W)
	Vr<--log(1-p)%*%t(W)
	Vbl<--log(bp)%*%t(W)
	Vbr<--log(1-bp)%*%t(W)
    V<-pmax(Vl,Vr)
    Vb<-pmax(Vbl,Vbr)
    nr<-nrow(V)
    np<-nrow(Vb)
    minV<-index.AW<-pval<-qval<-matrix(NA,nr,1)
    pV<-matrix(NA,nr,nW)
    minVb<-matrix(9999,nr,np/nr) 
    minVb[-rnum,]<-NA
    for (i in 1:nW) {
        m.Vb<-matrix(Vb[,i],nr,np/nr)
        pV[rnum,i]<-perm.p(as.matrix(V[rnum,i]),m.Vb[rnum,],tail="high")
        minVb<-pmin(minVb,emp(m.Vb[rnum,],tail="high"))
    }
    minV[rnum]<-apply(pV[rnum,],1,min)
	index.AW[rnum]<-apply(pV[rnum,],1,which.min)
	AW.weight<-W[index.AW,]
    pval[rnum]<-perm.p(minV[rnum],minVb[rnum,],tail="low")
    qval[rnum]<-p.adjust(pval[rnum],method="BH")
    row.names(minV)<-row.names(pval)<-row.names(qval)<-row.names(AW.weight)<-gene.names
    res<-list(stat=minV, pval=pval, FDR=qval, AW.weight=AW.weight)
    return(res)
}

cal.MCC<-function(dt1,dt2,l1,l2)
{
 l1<-unclass(factor(l1))
 l2<-unclass(factor(l2))
 K<-nlevels(l1)
 n1<-table(factor(l1,levels=unique(l1)))
 n2<-table(factor(l2,levels=unique(l2)))
 ind1<-diag(rep(1,length(n1)))[rep(1:nlevels(l1),n1),] 
 ind2<-diag(rep(1,length(n2)))[rep(1:nlevels(l2),n2),]
 xk.<-dt1%*%ind1%*%diag(1/n1)
 yk.<-dt2%*%ind2%*%diag(1/n2)
 x..<-rowMeans(xk.)
 y..<-rowMeans(yk.)
 sxk.yk.<-rowSums(xk.*yk.)
 num<-1/K*sxk.yk.-x..*y..
 sumx2<-dt1^2%*%ind1
 sumy2<-dt2^2%*%ind2
 vx<-1/K*rowSums((sumx2-xk.^2)%*%diag(1/(n1-1)))-x..^2
 vy<-1/K*rowSums((sumy2-yk.^2)%*%diag(1/(n2-1)))-y..^2
 den<-sqrt(vx*vy)
 r<-num/den
 return(r)
}

cal.minMCC<-function(dat,lbl)
{
 K<-length(dat)  
 if (K==2) min.MCC<-cal.MCC(dat[[1]],dat[[2]],factor(lbl[[1]]),factor(lbl[[2]]))
 else
 {
  allcomb<-combn(1:K,2) 
  pair.mcc<-NULL 
  for (i in 1:ncol(allcomb))
  {
   dt1<-dat[[allcomb[1,i]]] 
   dt2<-dat[[allcomb[2,i]]]
   l1<-lbl[[allcomb[1,i]]]
   l2<-lbl[[allcomb[2,i]]]
   pair.mcc<-cbind(pair.mcc,cal.MCC(dt1,dt2,factor(l1),factor(l2)))
  }
   row.names(pair.mcc)<-row.names(dat[[1]])
   min.MCC<-apply(pair.mcc,1,min) 
  }
  return(min.MCC)
}


get.minMCC<-function(dat,lbl,nperm)
{
  gene.names<-row.names(dat[[1]]) 
  rnum<-unlist(lapply(dat,function(x) which(apply(x,1,function(x) any(!is.na(x))))))
   perm.mcc.perm<-function(dat,lbl)
   {
     perm.d<-list()
      for (i in 1:length(dat))
    {
      perm.d[[i]]<-dat[[i]][,sample(1:ncol(dat[[i]]))]
     }
     perm.r<-cal.minMCC(perm.d,lbl)
    return(perm.r)
    }
   Ug<-cal.minMCC(dat,lbl)
  if (!is.null(nperm))
    {
     Ubg<-replicate(nperm,perm.mcc.perm(dat=dat,lbl=lbl))
    }
   else
    {
     stop("there're no parametric results for MCC statistic,you need to specify a number to nperm")   
    }
     pval<-qval<-matrix(NA,nrow(dat[[1]]),1)
     pval[rnum]<-perm.p(Ug[rnum],Ubg[rnum],tail="high")
     qval[rnum]<-p.adjust(pval[rnum],method="BH")
     names(Ug)<-row.names(pval)<-row.names(qval)<-gene.names
     res<-list(stat=as.matrix(Ug),pval=as.matrix(pval),FDR=as.matrix(qval))
     attr(res,"nstudy")<-length(dat)
     attr(res,"meta.method")<-"minMCC"
     return(res)
}


perm.p<-function(stat,perm,tail) {
	G<-length(stat)
      B<-length(perm)/G
	if(tail=="low"){
		r = rank(c(stat, as.vector(perm)),ties.method="max")[1:G]
    		r2 = rank(c(stat),ties.method="max")
    		r3 = r - r2
    		p = r3/(B*G)
	}
	if(tail=="high"){
		r = rank(c(stat, as.vector(perm)),ties.method="min")[1:G]
    		r2 = rank(stat,ties.method="max")
    		r3 = r - r2
    		p = 1-r3/(B*G)
	}
	if(tail=="abs"){
		r = rank(c(abs(stat), abs(as.vector(perm))),ties.method="min")[1:G]
    		r2 = rank(c(abs(stat)),ties.method="max")
    		r3 = r - r2
    		p = 1-r3/(B*G)
	}
	p[p==0]<-1e-20
	p[p==1]<-1-1e-10
    	return(p)
}



emp<-function(mat,tail) {
	B<-ncol(mat)
    	G<-nrow(mat)
    	if (tail=="low"){
		s<-matrix(rank(mat,ties.method="max"),G,B)
		p<-s/G/B}
    	if (tail=="high"){
	 	s<-matrix(rank(mat,ties.method="min"),G,B)
		 p<-1-(s-1)/G/B}
   	 if (tail=="abs"){
		s<-matrix(rank(abs(mat),ties.method="min"),G,B)
		p<-1-(s-1)/G/B}

    p[p==0]<-1e-20
    p[p==1]<-1-1e-10
    return(p)
}


CDF.rop<-function(z,r,n){
	require(combinat)
	if (r>=.5*(n+1)){
		pval<-ifelse(z>=0.5&z<1,
			1-sum(sapply((n-r+1):(r-1),function(y)sum(sapply((n-r+1):(n-y),function(x)dmnom(c(y,x,n-y-x),n,c(1-z,1-z,2*z-1))))))				
			,2*(1-pbinom(r-1,n,z)))}
	else{
		pval<-ifelse(z>=0&z<=0.5,
			1-sum(sapply((0):(r-1),function(y)sum(sapply(0:(r-1),function(x)dmnom(c(y,x,n-y-x),n,c(z,z,1-2*z)))))),				
			1)}		
	return(pval)
	}


get.SR<-function(p,bp=NULL){
	k<-ncol(p)
	rnum<-which(apply(p,1,function(x) !any(is.na(x))))
	pval<-Ug<-rep(NA,nrow(p))
	Ug[rnum]<-rowSums(apply(p,2,rank))[rnum]
	if(!is.null(bp)){
		nperm<-nrow(bp)/nrow(p)
		Ubg<-matrix(NA,nrow(p),nperm)
		for(j in 1:nperm){
			Ubg[rnum,j]<-rowSums(apply(bp[((j-1)*nrow(p)+1):(j*nrow(p)),],2,rank))[rnum]
		}
		pval[rnum]<-perm.p(Ug[rnum],Ubg[rnum,],tail="low")
		qval<-p.adjust(pval,method="BH")
	}else{
		nperm=500
		bp<-matrix(runif(500*nrow(p)*k),500*nrow(p),k)
		Ubg<-matrix(NA,nrow(p),nperm)
		for(j in 1:nperm){
			Ubg[rnum,j]<-rowSums(apply(bp[((j-1)*nrow(p)+1):(j*nrow(p)),],2,rank))[rnum]
		}
		pval[rnum]<-perm.p(Ug[rnum],Ubg[rnum,],tail="low")
		qval<-p.adjust(pval,method="BH")}
	res<-list(stat=Ug, pval=pval, FDR=qval)	
	names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
    return(res)
}


get.PR<-function(p,bp=NULL){
	rowProds<-function (x, ...){
    		s <- (x == 0)
    		s <- rowSums(s,na.rm=T)
   		 ok <- (s == 0)
    		rm(s)
   		 x <- x[ok, , drop = FALSE]
    		y <- vector(mode(x), nrow(x))
    		s <- (x < 0)
    		s <- rowSums(s,na.rm=T)
    		s <- (s%%2)
    		s <- c(+1, -1)[s + 1]
    		x <- abs(x)
    		x <- log(x)
    		x <- rowSums(x, ...)
    		x <- exp(x)
    		x <- s * x
    		y[ok] <- x
   		 rm(ok, s, x)
   		 y
	}
	k<-ncol(p)
	pval<-Ug<-rep(NA,nrow(p))
	rnum<-which(apply(p,1,function(x) !any(is.na(x))))
	Ug[rnum]<-rowProds(apply(p,2,rank),na.rm=T)[rnum]
	if(!is.null(bp)){
		nperm<-nrow(bp)/nrow(p)
		Ubg<-matrix(NA,nrow(p),nperm)
		for(j in 1:nperm){
			Ubg[rnum,j]<-rowProds(apply(bp[((j-1)*nrow(p)+1):(j*nrow(p)),],2,rank),na.rm=T)[rnum]
		}
		pval[rnum]<-perm.p(Ug[rnum],Ubg[rnum,],tail="low")
		qval<-p.adjust(pval,method="BH")
	}else{
		nperm=500
		bp<-matrix(runif(500*nrow(p)*k),500*nrow(p),k)
		Ubg<-matrix(NA,nrow(p),nperm)
		for(j in 1:nperm){
			Ubg[rnum,j]<-rowProds(apply(bp[((j-1)*nrow(p)+1):(j*nrow(p)),],2,rank),na.rm=T)[rnum]
		}
		pval[rnum]<-perm.p(Ug[rnum],Ubg[rnum,],tail="low")
		qval<-p.adjust(pval,method="BH")}
	res<-list(stat=Ug, pval=pval, FDR=qval)	
	names(res$stat)<-names(res$pval)<-names(res$FDR)<-rownames(p)
    return(res)
}


cal.ES<-function(y,l,paired=FALSE){
	l<-unclass(factor(l))
	n<-table(factor(l))
	if(paired){
		if (n[1]!=n[2]) {
            stop("The study is not paired design")
        }
		DM<-y[,l==2]
		CM<-y[,l==1]
		ydiff<-DM-CM
		den<-sqrt(1/(n[1]-1)*(rowSums(ydiff^2))-1/(n[1]^2-n[1])*(rowSums(ydiff))^2)
		t<-rowMeans(ydiff)/(den/sqrt(n[1]))
		rnum<-n[1]*rowSums(DM*CM)-rowSums(DM)*rowSums(CM)
		rden<-sqrt((n[1]*rowSums(DM^2)-(rowSums(DM))^2)*(n[1]*rowSums(CM^2)-(rowSums(CM))^2))
		r<-rnum/rden
		d<-t*sqrt(2*(1-r)/n[1])
		m<-n[1]-1
		cm = gamma(m/2)/(sqrt(m/2) * gamma((m - 1)/2))
		dprime=cm*d
		vard=(2*(1-r)/n[1])*((n[1]-1)/(n[1]-3))*(1+n[1]*dprime^2/(2*(1-r)))-dprime^2/cm^2
		vardprime=cm^2*vard
	
	}else{
		ind<-diag(rep(1,length(n)))[l,]
		ym<-y%*%ind%*%diag(1/n) 
		ntilde<-1/sum(1/n)
		m=sum(n)-2
		cm = gamma(m/2)/(sqrt(m/2) * gamma((m - 1)/2))
		s<-sqrt((1/(sum(n)-2)*((y^2%*%ind)%*%diag(1/(n-1))-ym^2%*%diag(n/(n-1)))%*%(n-1)))
		d<-(ym[,2]-ym[,1])/s
		dprime=d-3*d/(4*(sum(n)-2)-1)
		terme1=1/ntilde
		vard = terme1 + d^2 * (terme1 * ntilde - 1/cm^2)
		vardprime=sum(1/n)+dprime^2/(2*sum(n))			
	}
	result = cbind( dprime, vardprime)
    colnames(result) = c( "dprime", "vardprime")
	rownames(result)<-rownames(y)
    result
}

get.ES<-function(x,paired){   
	K<-length(x) 
	ES.m<-Var.m<- N<-n<-NULL  
	for (k in 1:K){
		y<-x[[k]][[1]]
		l<-x[[k]][[2]]
		temp<-cal.ES(y,l,paired[k])
		ES.m<-cbind(ES.m,temp[,"dprime"])
		Var.m<-cbind(Var.m,temp[,"vardprime"])
	     N<-c(N,length(l))
	     n<-c(n,table(l))
	}
	rownames(ES.m)<-rownames(y)
	rownames(Var.m)<-rownames(y)
	colnames(ES.m)<-paste("study",1:K,sep="")
	colnames(Var.m)<-paste("study",1:K,sep="")
	res<-list(ES=ES.m,Var=Var.m)
      attr(res,"nperstudy")<-N
      attr(res,"nperlabelperstudy")<-n 
	return(res)	
}

ind.cal.ES<-function(x,paired,nperm=NULL,miss.tol=0.3){
	x<-MetaDE.impute(x,y=miss.tol)
	K<-length(x)
	res<-get.ES(x,paired=paired)
	if(!is.null(nperm)){
		perm.ES<-perm.Var<-NULL
		for(i in 1:nperm){
			for(k in 1:K){
				x[[k]][[2]]<-perm.lab(x[[k]][[2]],paired[k])
			}
			tempRes<-get.ES(x,paired=paired)
			perm.ES<-rbind(perm.ES,tempRes$ES)
			perm.Var<-rbind(perm.Var,tempRes$Var)
		}
	}else{
		perm.ES<-perm.Var<-NULL
	}
	if(is.null(names(x))){colnames(res$ES)<-colnames(res$Var)<-paste("dataset",1:K,sep="")
	}else{colnames(res$ES)<-colnames(res$Var)<-names(x)}
      result<-list(ES=res$ES,Var=res$Var,perm.ES=perm.ES,perm.Var=perm.Var)
      attr(result,"nperstudy")<-attr(res,"nperstudy")
      attr(result,"nperlabelperstudy")<-attr(res,"nperlabelperstudy") 
	return(result)
}


get.Q<-function(em,vm){
	wt <- 1/vm
    temp1 <- wt * em
    mu.hat <- rowSums(temp1)/rowSums(wt)
    Q <- rowSums(wt * (em - mu.hat)^2)
	return(Q)
}

get.tau2<-function(Q,vm,k){
	wt<-1/vm
	s1 <- rowSums(wt)
    s2 <- rowSums(wt^2)
    temp<- (Q - (k - 1))/(s1 - (s2/s1))
    tau2<-pmax(temp,0)
	return(tau2)
}


get.FEM<-function(em,vm,pe=NULL,pv=NULL){
	wt<-1/vm
	mu.hat<-rowSums(wt*em)/rowSums(wt)
	mu.var<-1/rowSums(wt)
	z.score<-mu.hat/sqrt(mu.var)
	if(!is.null(pe)&!is.null(pv)){
		rnum<-which(apply(em,1,function(x) !any(is.na(x))))
		Z0<-matrix(get.REM2(pe,pv)$zval,nrow(em),nrow(pe)/nrow(em))
		z.p<-rep(NA,nrow(em))
		z.p[rnum]<-perm.p(z.score[rnum],Z0[rnum,],"abs")
	}else{
		z.p<-2*(1-pnorm(abs(z.score)))
	}
	qval<-p.adjust(z.p,method="BH")
	res<-list(mu.hat=mu.hat,mu.var=mu.var,zval=z.score,pval=z.p,FDR=qval)
	return(res)
}

get.FEM2<-function(em,vm){
	wt<-1/vm
	mu.hat<-rowSums(wt*em)/rowSums(wt)
	mu.var<-1/rowSums(wt)
	z.score<-mu.hat/sqrt(mu.var)
	z.p<-2*(1-pnorm(abs(z.score)))
	qval<-p.adjust(z.p,method="BH")
	res<-list(mu.hat=mu.hat,mu.var=mu.var,zval=z.score,pval=z.p,FDR=qval)
	return(res)
}

get.REM2<-function(em,vm){
	k<-ncol(em)
	Q.val<-get.Q(em,vm)
	tau2<-get.tau2(Q.val,vm,k)
	temp.wt<-1/(vm+tau2)	
	mu.hat<-rowSums(temp.wt*em)/rowSums(temp.wt)
	mu.var<-1/rowSums(temp.wt)
	Qpval <- pchisq(Q.val, df = k - 1, lower.tail = FALSE)
	z.score<-mu.hat/sqrt(mu.var)
	z.p<-2*(1-pnorm(abs(z.score)))
	qval<-p.adjust(z.p,method="BH")
	res<-list(mu.hat=mu.hat,mu.var=mu.var,Qval=Q.val,Qpval=Qpval,tau2=tau2,zval=z.score,pval=z.p,FDR=qval)
	return(res)
}

get.REM<-function(em,vm,pe=NULL,pv=NULL){
	k<-ncol(em)
	Q.val<-get.Q(em,vm)
	tau2<-get.tau2(Q.val,vm,k)
	temp.wt<-1/(vm+tau2)	
	mu.hat<-rowSums(temp.wt*em)/rowSums(temp.wt)
	mu.var<-1/rowSums(temp.wt)
	Qpval <- pchisq(Q.val, df = k - 1, lower.tail = FALSE)
	z.score<-get.REM2(em,vm)$zval
	if(!is.null(pe)&!is.null(pv)){
		rnum<-which(apply(em,1,function(x) !any(is.na(x))))
		Z0<-matrix(get.REM2(pe,pv)$zval,nrow(em),nrow(pe)/nrow(em))
		z.p<-rep(NA,nrow(em))
		z.p[rnum]<-perm.p(z.score[rnum],Z0[rnum,],"abs")
	}else{
		z.p<-2*(1-pnorm(abs(z.score)))
	}
	qval<-p.adjust(z.p,method="BH")
	res<-list(mu.hat=mu.hat,mu.var=mu.var,Qval=Q.val,Qpval=Qpval,tau2=tau2,zval=z.score,pval=z.p,FDR=qval)
	return(res)
}

get.RP<-function(dat,lbl, nperm = 100, logged = TRUE) 
{
  gene.names<-row.names(dat[[1]]) 
  num.perm<-nperm
  num.ori<-length(lbl) 
  num.gene<-nrow(dat[[1]]) 
  nk<-unlist(lapply(lbl,function(x) length(x))) 
  origin<-rep(1:num.ori,nk)
    
  data<-cl<-NULL
  for (k in 1:num.ori)
  {
   data<-cbind(data,dat[[k]])
   cl<-c(cl,lbl[[k]])
  }
  
    total.sam = length(origin)
    total.sam1<-ncol(data)
    if (total.sam != total.sam1) 
        stop("the lbl number does not match the dat columns")

    data.pre = OriginxyCall(data=data, cl=cl, origin=origin)
    y = data.pre$data2
   
    data = as.matrix(data)
    mode(data) = "numeric"
    NA.genes <- NULL
    if (any(is.na(data))) {
        NA.genes <- unique(ceiling(which(is.na(t(data)))/ncol(data)))
        cat("Warning: There are", length(NA.genes), "genes with at least one missing value.", 
            "\n", "\n")
        if (na.rm) 
            data[NA.genes, ] <- NaReplace2(data[NA.genes, ],origin)
        if (!na.rm) 
            cat(" This value is not used to compute rank product.", 
                "\n", "\n")
    }

    if (!is.null(y)) {
        num.class = 2
        data1.all = data.pre$data1
        data2.all = data.pre$data2
        fold.change = NULL
        for (l in 1:num.ori) {
            data1 = as.matrix(data1.all[[l]])
            data2 = as.matrix(data2.all[[l]])
            data1.ave = apply(data1, 1, mean)
            data2.ave = apply(data2, 1, mean)
            if (logged) {
                fold.change = cbind(fold.change,(data1.ave-data2.ave))
            }
            else {
                fold.change = cbind(fold.change, (data1.ave/data2.ave))
            }
            rm(data1, data2, data1.ave, data2.ave)
        }
        ave.fold.change = apply(fold.change, 1, mean)
    }
    if (is.null(y)) {
        num.class = 1
        data1.all = data.pre$data1
        data2.all = data.pre$data2
        fold.change = NULL
        for (l in 1:num.ori) {
            data1 = as.matrix(data1.all[[l]])
            fold.change = cbind(fold.change, apply(data1, 1, 
                mean))
            rm(data1)
        }
        ave.fold.change = apply(fold.change, 1, mean)
    }
    RP.ori.out.upin2 = RankProd2(data1.all, data2.all, num.ori, 
        num.gene, logged, num.class, rev.sorting = FALSE)
    RP.ori.upin2 = RP.ori.out.upin2$RP       
    rank.ori.upin2 = rank(RP.ori.upin2)
    RP.ori.out.downin2 = RankProd2(data1.all, data2.all, num.ori, 
        num.gene, logged, num.class, rev.sorting = TRUE)
    RP.ori.downin2 = RP.ori.out.downin2$RP   
    rank.ori.downin2 = rank(RP.ori.downin2)
    RP.perm.upin2 <- matrix(NA, num.gene, num.perm)
    RP.perm.downin2 <- matrix(NA, num.gene, num.perm)
    cat("Starting ", num.perm, "permutations...", "\n")
    for (p in 1:num.perm) {
        new.data.temp = NewdataCom(data1.all, data2.all, num.ori,num.class)
        new.data1.all = new.data.temp$new.data1.all
        new.data2.all = new.data.temp$new.data2.all
        temp1 = RankProd2(new.data1.all, new.data2.all, num.ori,num.gene, logged, num.class, rev.sorting = FALSE)
        RP.perm.upin2[, p] = temp1$RP
        rm(temp1)
        temp2 = RankProd2(new.data1.all, new.data2.all, num.ori,num.gene, logged, num.class, rev.sorting = TRUE)
        RP.perm.downin2[, p] = temp2$RP
        rm(temp2)
    }
    pval.upin2<-perm.p(RP.ori.upin2,RP.perm.upin2,tail="low") 
    pval.downin2<-perm.p(RP.ori.downin2,RP.perm.downin2,tail="low") 
    qval.upin2<-p.adjust(pval.upin2,method="BH")
    qval.downin2<-p.adjust(pval.downin2,method="BH")
    names(RP.ori.upin2)<-names(pval.upin2)<-names(qval.upin2)<- names(RP.ori.downin2)<-names(pval.downin2)<-names(qval.downin2)<-gene.names

    res<-list(meta.stat.up=RP.ori.upin2,pval.up=pval.upin2,FDR.up=qval.upin2,
              meta.stat.down=RP.ori.downin2,pval.down=pval.downin2,FDR.down=qval.downin2,
    AveFC = ave.fold.change)
    return(res)
}

OriginxyCall<-function (data, cl, origin, sum = FALSE) 
{
    lev <- unique(cl)
    uni.cl <- length(lev)
    if (uni.cl > 2) 
        stop("There is something wrong with the classlabels")
    ori.lev <- unique(origin)
    uni.ori <- length(ori.lev)
    cat(" The data is from ", uni.ori, "different origins \n \n")
    if (min(ori.lev) != 1 | max(ori.lev) != uni.ori) {
        cat("Warning: origins labels are not numbered from 1 to ", 
            uni.ori, "\n", "\n")
    }
    if (uni.cl == 1) {
        if (sum) {
            cat("Rank Sum analysis for one-class case", "\n", 
                "\n")
        }
        else {
            cat("Rank Product analysis for one-class case", "\n", 
                "\n")
        }
        if (lev != 1) {
            cat("warning: Expected classlabel is 1, cl will hus be set to 1.", 
                "\n", "\n")
            cl = rep(1, length(cl))
        }
        data2 <- NULL
        data1 <- vector("list", uni.ori)
        for (c in 1:uni.ori) {
            data1[[c]] = data[, origin == ori.lev[[c]]]
        }
    }
    if (uni.cl == 2) {
        if (sum) {
            cat("Rank Sum analysis for two-class case", "\n", 
                "\n")
        }
        else {
            cat("Rank Product analysis for two-class case", "\n", 
                "\n")
        }
        if (min(lev) != 0 | max(lev) != 1) {
            cat("Warning: Expected classlabels are 0 and 1. cl will thus be set to 0 and 1.", 
                "\n", "\n")
            cl[which(cl == min(lev))] <- 0
            cl[which(cl == max(lev))] <- 1
        }
        data2 <- vector("list", uni.ori)
        data1 <- vector("list", uni.ori)
        for (c in 1:uni.ori) {
            index1 <- which(origin == ori.lev[[c]] & cl == 0)
            index2 <- which(origin == ori.lev[[c]] & cl == 1)
            if (length(index1) == 0 | length(index1) == 0) 
                stop("Error: data from different origins should contain data from both classs")
            data1[[c]] <- data[, index1]
            data2[[c]] <- data[, index2]
            rm(index1, index2)
        }
    }
    list(data1 = data1, data2 = data2)
}

RankProd2<-function (data1.all, data2.all, num.ori, num.gene, logged, num.class, 
    rev.sorting) 
{
    num.rank.all = 0
    rank.rep.all = t(t(1:num.gene))
    for (l in 1:num.ori) {
        data1 = data1.all[[l]]
        data2 = data2.all[[l]]
        data1 = as.matrix(data1)
        if (num.class == 2) {
            data2 = as.matrix(data2)
        }
        temp = RankComp(data1, data2, logged, num.class, rev.sorting)
        rank.rep.all = cbind(rank.rep.all, temp$rank)
        num.rank.all = num.rank.all + temp$num.rank
        rm(temp)
    }
    rank.all = rank.rep.all[, -1]
    rank.prod.temp = rank.all^(1/num.rank.all)
    rank.prod = apply(rank.prod.temp, 1, prod)
    rank.prod[num.rank.all == 0] = NA
    list(RP = rank.prod, rank.all = rank.all)
}
RankComp<-function (data1, data2, logged, num.class, rev.sorting) 
{
    num.gene = dim(data1)[1]
    if (num.class == 2) {
        if (rev.sorting) {
            data1.wk <- data2
            data2.wk <- data1
        }
        else {
            data1.wk <- data1
            data2.wk <- data2
        }
        k1 = dim(data1.wk)[2]
        k2 = dim(data2.wk)[2]
        num.rep = k1 * k2
        data.rep = matrix(NA, num.gene, num.rep)
        if (logged) {
            for (k in 1:k1) {
                temp = ((k - 1) * k2 + 1):(k * k2)
                data.rep[, temp] = data1.wk[, k] - data2.wk
            }
        }
        else {
            for (k in 1:k1) {
                temp = ((k - 1) * k2 + 1):(k * k2)
                data.rep[, temp] = data1.wk[, k]/data2.wk
            }
        }
        rank.rep = apply(data.rep, 2, rank)
    }
    if (num.class == 1) {
        data.rep = data1
        if (rev.sorting) {
            num.rep = dim(data1)[2]
            rank.temp = matrix(NA, num.gene, num.rep)
            for (r in 1:num.rep) {
                rank.temp[, r] = rank(data1[, r], na.last = FALSE)
            }
            rank.rep = (num.gene + 1) - rank.temp
        }
        else {
            rank.rep = apply(data1, 2, rank)
        }
    }
    rank.rep[is.na(data.rep)] = 1
    num.rank = apply(is.na(data.rep) == FALSE, 1, sum)
    list(rank = rank.rep, num.rank = num.rank)
}

NewdataCom <-function(data1.all,data2.all,num.ori,num.class)
{
  new.data1.all <- vector("list",num.ori)
  new.data2.all <- vector("list",num.ori)  
  for ( l in 1:num.ori ){
      
      data1 <- as.matrix(data1.all[[l]])
      if (num.class == 2) {data2 <- as.matrix(data2.all[[l]])}

      temp.data <- Newdata(data1,data2,num.class)
      new.data1.all[[l]] <- temp.data$new.data1
      new.data2.all[[l]] <- temp.data$new.data2
  }
  if(num.class == 1) {new.data2.all <- NULL} 
  list(new.data1.all = new.data1.all,new.data2.all = new.data2.all)
}

Newdata<-function(data1,data2,num.class)
{
  
  k1 <- dim(data1)[2]
  num.gene <- dim(data1)[1]
  new.data1 <- matrix(NA,num.gene,k1)
  
  for (k_1 in 1:k1)
  {
    temp <- sample(1:num.gene,num.gene,replace = FALSE, prob = NULL)
    new.data1[,k_1] <- data1[temp,k_1];
    rm(temp)
   }
   rm(k_1,k1)

  if (num.class == 2) {
     k2 <- dim(data2)[2]
     new.data2 <- matrix(NA,num.gene,k2) 
     for (k_2 in 1:k2) {   
         temp <- sample(1:num.gene,num.gene,replace = FALSE, prob = NULL)
         new.data2[,k_2] <- data2[temp,k_2];
     }
      rm(k_2,k2)
   }

   if (num.class == 1) { new.data2=NULL}
   list(new.data1 = new.data1,new.data2 = new.data2)
   
 }