# Package: bmem
# Type: Package
# Title: Mediation analysis with missing data using bootstrap
# Version: 1.3
# Date: 2011-01-04
# Author: Zhiyong Zhang and Lijuan Wang
# Maintainer: Zhiyong Zhang <zhiyongzhang@nd.edu>
# Depends: R (>= 1.7), sem, Amelia, MASS
# Description: Four methods for mediation analysis with missing data: Listwise deletion, Pairwise deletion, Multiple imputation, and Two Stage Maximum Likelihood algorithm. For MI and TS-ML, auxiliary variables can be included. Bootstrap confidence intervals for mediation effects are obtained. The robust method is also implemented for TS-ML.


bmem.moments<-function(x, type=0){
	#0: listwise deletion
	#1: pairwise deletion
	x<-cbind(1, x)
	n<-nrow(x)
	p<-ncol(x)
	res<-matrix(NA, p, p)	
	x<-as.matrix(x)
	if (type ==0){
		x<-na.omit(x)
		res<-crossprod(x,x)/n
	}else{
		for (i in 1:p) res[i,i]<-mean( (x[complete.cases(x[,i]),i])^2 )
		for (i in 1:(p-1)){
			for (j in (i+1):p){
				y<-x[,c(i,j)]
				y<-na.omit(y)
				res[i,j]<-res[j,i]<-mean(y[,1]*y[,2])
			}
		}
	}
	rownames(res)<-colnames(res)<-c('(i)', (colnames(x))[2:p])
	res[1,1]<-1
	res
}

bmem.list.cov<-function(x, moment=FALSE){
	## x: data
	## moment=FALSE: covariance analysis only
	listdata<-na.omit(x)	
	if (moment){ 
		temp.cov<-bmem.moments(listdata, type=0)
	}else{
		temp.cov<-cov(listdata)
	}
	temp.cov
}

bmem.pair.cov<-function(x, moment=FALSE){
	## x: data
	## moment=FALSE: covariance analysis only	
	if (moment){
		temp.cov<-bmem.moments(x, type=1)
	}else{
		temp.cov<-cov(x, use='pairwise.complete.obs')
	}
	temp.cov
}

bmem.sem<-function(x, ram, N, indirect, moment=FALSE, ...){
	if (moment){
	    sem.res<-sem(ram, x, N, raw=TRUE, fixed.x='(i)', ...)
		est<-sem.res$coeff
	}else{
	    sem.res<-sem(ram, x, N, ...)
		est<-sem.res$coeff
	}
	
	if (!missing(indirect)){
	n<-length(indirect)
	est.indirect<-rep(0,n)
	for (i in 1:n){
		temp<-indirect[i]
		temp<-gsub(' +','',temp)
		temp<-gsub('-','+', temp, fixed=TRUE)
		temp<-unlist(strsplit(temp, '+', fixed=TRUE))
		m<-length(temp)
		temp.est<-0
		par<-NULL
		for (j in 1:m){
			temp1<-unlist(strsplit(temp[j], '*', fixed=TRUE))
			par<-c(par, temp1)
		}
		ind.exp<-parse(text=indirect[i])
		par.list<-as.list(est[par])
		est.indirect[i]<-eval(ind.exp, par.list)
		}
	names(est.indirect)<-indirect
	}else{
		est.indirect<-NULL
	}
	model.fit<-NA
	model.summary<-summary(sem.res)
	model.fit<-c(model.summary$chisq,model.summary$GFI, model.summary$AGFI, model.summary$RMSEA[1], model.summary$NFI, model.summary$NNFI, model.summary$CFI, model.summary$BIC, model.summary$SRMR)
	
	list(est=c(est, est.indirect), model.fit=model.fit)
	}

bmem.sobel.ind<-function(sem.object, ind){
	est<-sem.object$coeff
	temp<-gsub(' +','',ind)	
        temp<-gsub('-','+', temp, fixed=TRUE)
	temp<-unlist(strsplit(temp, '+', fixed=TRUE))
	m<-length(temp)
	temp.est<-0
	par<-NULL
	for (j in 1:m){
		temp1<-unlist(strsplit(temp[j], '*', fixed=TRUE))
		par<-c(par, temp1)
		}
	ind.exp<-parse(text=ind)
	par.list<-as.list(est[par])	
	est.indirect<-eval(ind.exp, par.list)
	ind.deriv<-deriv(ind.exp, par)
	first.deriv<-eval(ind.deriv, par.list)
	first.deriv<-attributes(first.deriv)$gradient
	
	if (is.null(sem.object$cov)) {
		var.par<-sem.object$vcov[par,par]
		}else{
			var.par<-sem.object$cov[par,par]
		}
	var.ind<-first.deriv%*%var.par%*%t(first.deriv)
	s.e.ind<-sqrt(var.ind)
	
	res<-matrix(c(est.indirect, s.e.ind, est.indirect/s.e.ind), 1, 3)
	colnames(res)<-c('Estimate', 'S.E.', 'z-score')
	rownames(res)<-ind
	res
	}

bmem.sobel<-function(x, ram, indirect, moment=FALSE, ...){
	N<-nrow(x)
	temp.cov<-bmem.list.cov(x, moment)
	if (moment){
	    sem.object<-sem(ram, temp.cov, N, raw=TRUE, fixed.x='(i)', ...)
	}else{
	    sem.object<-sem(ram, temp.cov, N, ...)
	}
	if (!missing(indirect)){
		n<-length(indirect)
		est.indirect<-NULL
		for (i in 1:n){
			temp<-indirect[i]
			res<-bmem.sobel.ind(sem.object, temp)
			est.indirect<-rbind(est.indirect, res)
		}
		sem.est<-summary(sem.object)$coeff[,1:3]
		colnames(sem.est)<-colnames(est.indirect)
	}else{
		est.indirect<-NULL
		sem.est<-summary(sem.object)$coeff[,1:3]
	}
		
	#print(summary(sem.object))
	all.res<-rbind(sem.est, est.indirect)
	p.value<-2*(1-pnorm(abs(all.res[,3])))
	all.res<-cbind(all.res, p.value)
			
	##cat('\n\nIndirect or mediation effects\n')
	print(all.res)
	invisible(list(estimates=all.res, sem=sem.object))
	}


bmem.list<-function(x, ram, indirect, moment=FALSE, ...){
	N<-dim(x)[1]
	temp.cov<-bmem.list.cov(x, moment)
	bmem.sem(temp.cov, ram, N, indirect, moment, ...)
	}

bmem.pair<-function(x, ram, indirect, moment=FALSE, ...){
	N<-dim(x)[1]
	temp.cov<-bmem.pair.cov(x, moment)
	bmem.sem(temp.cov, ram, N, indirect, moment, ...)
	}

bmem.mi.cov<-function(x, m=10, moment=FALSE){
	mi.data<-amelia(x,m,0)$imputations
	mi.cov<-list()
	for (j in 1:m){
		temp.cov<-bmem.pair.cov(mi.data[[j]], moment)
		mi.cov[[j]]<-temp.cov
		}
	mi.cov	
	}

bmem.v<-function(x, v, moment=FALSE){
	if (moment){
			v<-v+1
			v<-c(1,v)			
		}
		x[v,v]
	}

bmem.mi<-function(x, ram, indirect, v, m=10, moment=FALSE, ...){
	N<-dim(x)[1]
	if (missing(v)) v<-1:(ncol(x))
	temp.cov<-bmem.mi.cov(x,m,moment)
	mi.est<-NULL
	mi.fit<-NULL
	for (j in 1:m){
		v.cov<-bmem.v(temp.cov[[j]], v, moment)
		temp<-bmem.sem(v.cov, ram, N, indirect, moment, ...)
		mi.est<-rbind(mi.est, temp$est)
		mi.fit<-rbind(mi.fit, temp$model.fit)
		}
	mi.est<-apply(mi.est,2,mean)
	mi.fit<-apply(mi.fit,2,mean)
	list(est=mi.est, model.fit=mi.fit)
	}

bmem.pattern<-function(x){
	if (!is.matrix(x)) x<-data.matrix(x)
	n<-dim(x)[1]
	p<-dim(x)[2]
	misorder<-rep(0,n)
	for (i in 1:n){
		misorderj<-0
		for (j in 1:p){
			if (is.na(x[i,j])) misorderj<-misorderj+2^(j-1)
		}
		misorder[i]<-misorderj
	}
	temp<-order(misorder)
	x<-x[temp,]
	misn<-misorder[temp]
	
	mi<-0; nmi<-0;oi<-0; noi<-0;
	for (j in 1:p){
		if (is.na(x[1,j])){
			mi<-c(mi,j)  ## recording the missing variable subscript in the first case
			nmi<-nmi+1   ## number of missing values in the first case
		}else{
			oi<-c(oi,j)
			noi<-noi+1
		}
	}
	oi<-oi[2:(noi+1)]
	if (nmi==0){
		misinfo_0 = c(noi, oi)
	}else{
		mi<-mi[2:(nmi+1)]
		misinfo_0<-c(noi,oi,mi) ##recording the # observed variables, locations;
	}	
	patnobs <- 0 ## number of observed cases in a pattern
	totpat<-1; ncount<-1;
	t1<-misn[1]
	for (i in 2:n){
		if (misn[i]==t1){
			ncount<-ncount+1
		}else{
			patnobs<-c(patnobs,ncount)
			t1<-misn[i]
			ncount<-1
			totpat<-totpat+1
			mi<-0; nmi<-0;oi<-0; noi<-0;
			for (j in 1:p){
				if (is.na(x[i,j])){
					mi<-c(mi,j)
					nmi<-nmi+1
				}else{
					oi<-c(oi,j)
					noi<-noi+1
				}
			}
			oi<-oi[2:(noi+1)]
			mi<-mi[2:(nmi+1)]
			misinfo_0 <- rbind(misinfo_0, c(noi,oi,mi))
		}
	}
	patnobs<-c(patnobs, ncount)
	patnobs<-patnobs[2:(totpat+1)]
	if (is.vector(misinfo_0)) {misinfo<-c(patnobs, misinfo_0)}else{misinfo<-cbind(patnobs, misinfo_0)}
	if (!is.matrix(misinfo)){misinfo<-matrix(misinfo, nrow=1)}
	colnames(misinfo)<-NULL
	rownames(misinfo)<-NULL
	list(misinfo=misinfo, x=x)
}

bmem.ssq<-function(x){
	sum(x^2)
}

bmem.em.cov<-function(xmis, moment=FALSE, max_it=500){
	x<-xmis[[2]]
	misinfo<-xmis[[1]]
	ep <- 1e-12  ## precision
	n<-dim(x)[1]
	p<-dim(x)[2]
	mu0<-apply(x,2,mean, na.rm=TRUE) ## starting values for mu
	sig0<-cov(x,use="complete.obs") ## starting values for sigma
	n_it<-0;        
	err<-0
	dt<-1;
	while (dt>ep && n_it <= max_it){
		sumx<-rep(0,p); sumxx<-array(0,dim=c(p,p)); sumw1<-0; sumw2<-0;
		npat<-dim(misinfo)[1]  ## number of missing data patterns
		p1<-misinfo[1,2]       ## number of observed variables in pattern 1
		n1<-misinfo[1,1]       ## number of cases in pattern 1
		if (p1==p){            ## complete data
			sigin <- solve(sig0)  ## matrix inverse
			for (i in 1:n1){
				xi<-x[i,]
				xi0<-xi-mu0
				sumw1<-sumw1+1;
				xxi0<-xi0%*%t(xi0)
				sumx<-sumx+xi
				sumxx<-sumxx+xxi0
				sumw2<-sumw2+1
			} ## end for
		}else{ ## end p1==p
## with missing data			
			o1<-misinfo[1,3:(2+p1)]
			m1<-misinfo[1,(2+p1+1):(p+2)]
			mu_o<-mu0[o1]; mu_m<-mu0[m1]
			sig_oo<-sig0[o1,o1]; sig_om<-sig0[o1,m1];
			if (p1==1) {sig_mo<-sig_om}else{sig_mo<-t(sig_om)}
			sig_mm<-sig0[m1,m1];
			sigin_oo<-solve(sig_oo)
			beta_mo<-sig_mo%*%sigin_oo
			
			delt <- array(0, dim=c(p,p))
			delt[m1,m1]<-sig_mm - beta_mo%*%sig_om
			for (i in 1:n1){
				xi<-x[i,]
				xi_o<-xi[o1]
				xi0_o<-xi_o-mu_o
				stdxi_o<-sigin_oo%*%xi0_o
				sumw1<-sumw1+1
				xm1<-mu_m+sig_mo%*%stdxi_o
				xi[m1]<-xm1
				xi0<-xi-mu0
				xxi0<-xi0%*%t(xi0)
				sumx<-sumx+xi
				sumxx<-sumxx+xxi0+delt
				sumw2<-sumw2+1
			} ##end for 1:n1  
		}## end of (p1=p)
## start from pattern 2	
		if (npat>1){
			snj<-n1	
			for (j in 2:npat){
				nj<-misinfo[j,1]; pj<-misinfo[j,2];
				oj<-misinfo[j, 3:(2+pj)]; mj<-misinfo[j, (2+pj+1):(p+2)];
				mu_o<-mu0[oj]; mu_m<-mu0[mj];
				sig_oo<-sig0[oj,oj]; sig_om<-sig0[oj,mj];
				if (pj==1) {sig_mo<-sig_om}else{sig_mo<-t(sig_om)}  
				sig_mm<-sig0[mj,mj];
				sigin_oo<-solve(sig_oo)
				beta_mo<-sig_mo%*%sigin_oo
				delt <- array(0, dim=c(p,p))
				delt[mj,mj]<-sig_mm - beta_mo%*%sig_om
				
				for (i in ((snj+1):(snj+nj))){
					xi<-x[i,]
					xi_o<-xi[oj]
					xi0_o<-xi_o - mu_o
					stdxi_o<-sigin_oo%*%xi0_o					
					sumw1<-sumw1+1
					xmj<-mu_m+sig_mo%*%stdxi_o
					xi[mj]<-xmj
					xi0<-xi-mu0
					xxi0<-xi0%*%t(xi0)
					sumx<-sumx+1*xi
					sumxx<-sumxx+xxi0+delt
					sumw2<-sumw2+1
				}
				snj<-snj+nj
			} ## for (j in 2:npat)
		}
		mu1<-sumx/sumw1
		sig1<-sumxx/n
		dt<-(bmem.ssq(mu1-mu0)+bmem.ssq(sig1-sig0))/(bmem.ssq(mu0)+bmem.ssq(sig0));
		mu0<-mu1;
		sig0<-sig1;
		n_it<-n_it+1;
	} ## end while
	if (n_it == max_it) warnings('Maximum iteration for EM algorithm is exceeded!')
	temp.cov<-sig1
	rownames(temp.cov)<-colnames(sig1)
	m<-matrix(mu1, p, 1)
	if (moment){
		p1<-p+1
		temp.cov<-matrix(NA, p+1,p+1)
		temp.cov[1,1]<-1
		temp.cov[1,2:p1]<-mu1
		temp.cov[2:p1, 1]<-mu1
		temp.cov[2:p1, 2:p1]<-sig1+m%x%t(m)
		rownames(temp.cov)<-colnames(temp.cov)<-c('(i)', colnames(sig1))
	}
	temp.cov
}

bmem.em.rcov<-function(xmis, varphi=.1, moment=FALSE, max_it=1000, st='i'){
## x: data set
## misinfo: missing data pattern
## varphi: 
	x<-xmis[[2]]
	misinfo<-xmis[[1]]
	ep <- 1e-6  ## precision
	n<-dim(x)[1]
	p<-dim(x)[2]
	
	mu0<-rep(0,p)
	sig0<-diag(p)
	if (st=='mcd'){
		y<-na.omit(x)
		ny<-nrow(y)
		par.st<-cov.rob(y, method='mcd')
		mu0<-par.st$center
		sig0<-par.st$cov
	}

	n_it<-0;        
	dt<-1;
	if (varphi==0){
		ck<-10e+10
		cbeta<-1
	}else{
		prob<-1-varphi ## chi-square p-value
		chip<-qchisq(prob, p)
		ck<-sqrt(chip)
		cbeta<-( p*pchisq(chip, p+2) + chip*(1-prob) )/p
	}
	while (dt>ep && n_it <= max_it){
		sumx<-rep(0,p); sumxx<-array(0,dim=c(p,p)); sumw1<-0; sumw2<-0;
		npat<-dim(misinfo)[1]  ## number of missing data patterns
		p1<-misinfo[1,2]       ## number of observed variables in pattern 1
		n1<-misinfo[1,1]       ## number of cases in pattern 1
		if (p1==p){            ## complete data
			sigin <- solve(sig0)  ## matrix inverse
			for (i in 1:n1){
				xi<-x[i,]
				xi0<-xi-mu0
				di2<-xi0%*%sigin%*%xi0
				di<-sqrt(di2)
## Huber weight functions
				if (di<=ck){
					wi1<-1
					wi2<-1/cbeta
				}else{
					wi1<-ck/di
					wi2<-wi1*wi1/cbeta
				} ## end Huber weight
				sumw1<-sumw1+wi1;
				xxi0<-xi0%*%t(xi0)
				sumx<-sumx+wi1*xi
				sumxx<-sumxx+c(wi2)*xxi0
				sumw2<-sumw2+wi2
			} ## end for
		}else{ ## end p1==p
## with missing data
			if (varphi==0){
			   	ck1<-1e+10
			   	cbeta1<-1
			}else{ 
				chip1<-qchisq(prob, p1)
				ck1<-sqrt(chip1)
				cbeta1<-( p1*pchisq(chip1,p1+2) + chip1*(1-prob) )/p1
			}
			o1<-misinfo[1,3:(2+p1)]
			m1<-misinfo[1,(2+p1+1):(p+2)]
			mu_o<-mu0[o1]; mu_m<-mu0[m1]
			sig_oo<-sig0[o1,o1]; sig_om<-sig0[o1,m1];
			if (p1==1) {sig_mo<-sig_om}else{sig_mo<-t(sig_om)} 
			sig_mm<-sig0[m1,m1];
			sigin_oo<-solve(sig_oo)
			beta_mo<-sig_mo%*%sigin_oo
			
			delt <- array(0, dim=c(p,p))
			delt[m1,m1]<-sig_mm - beta_mo%*%sig_om
			for (i in 1:n1){
				xi<-x[i,]
				xi_o<-xi[o1]
				xi0_o<-xi_o-mu_o
				stdxi_o<-sigin_oo%*%xi0_o
				di2<-xi0_o%*%stdxi_o
				di<-sqrt(di2)
				if (di<=ck1){ ##Huber weight
					wi1<-1
					wi2<-1/cbeta1
				}else{
					wi1<-ck1/di
					wi2<-wi1*wi1/cbeta1
				}
				sumw1<-sumw1+wi1
				xm1<-mu_m+sig_mo%*%stdxi_o
				xi[m1]<-xm1
				xi0<-xi-mu0
				xxi0<-xi0%*%t(xi0)
				sumx<-sumx+wi1*xi
				sumxx<-sumxx+c(wi2)*xxi0+delt
				sumw2<-sumw2+wi2
			} ##end for 1:n1  
		}## end of (p1=p)
## start from pattern 2	
		if (npat>1){
			snj<-n1	
			for (j in 2:npat){
				nj<-misinfo[j,1]; pj<-misinfo[j,2];
				oj<-misinfo[j, 3:(2+pj)]; mj<-misinfo[j, (2+pj+1):(p+2)];
				mu_o<-mu0[oj]; mu_m<-mu0[mj];
				sig_oo<-sig0[oj,oj]; sig_om<-sig0[oj,mj];
				if (pj==1) {sig_mo<-sig_om}else{sig_mo<-t(sig_om)} 
				sig_mm<-sig0[mj,mj];
				sigin_oo<-solve(sig_oo)
				beta_mo<-sig_mo%*%sigin_oo
				delt <- array(0, dim=c(p,p))
				delt[mj,mj]<-sig_mm - beta_mo%*%sig_om
				if (varphi==0){
					ckj<-10e+10
					cbetaj<-1
				}else{
					chipj<-qchisq(prob,pj)
					ckj<-sqrt(chipj)
					cbetaj<- ( pj*pchisq(chipj, pj+2) + chipj*(1-prob) )/pj
				}
				for (i in ((snj+1):(snj+nj))){
					xi<-x[i,]
					xi_o<-xi[oj]
					xi0_o<-xi_o - mu_o
					stdxi_o<-sigin_oo%*%xi0_o
					di2<-xi0_o%*%stdxi_o
					di<-sqrt(di2)
					if (di<=ckj){ ##Huber weight
						wi1<-1
						wi2<-1/cbetaj
					}else{
						wi1<-ckj/di
						wi2<-wi1*wi1/cbetaj
					}
					sumw1<-sumw1+wi1
					xmj<-mu_m+sig_mo%*%stdxi_o
					xi[mj]<-xmj
					xi0<-xi-mu0
					xxi0<-xi0%*%t(xi0)
					sumx<-sumx+wi1*xi
					sumxx<-sumxx+c(wi2)*xxi0+delt
					sumw2<-sumw2+wi2
				}
				snj<-snj+nj
			} ## for (j in 2:npat)
		}
		mu1<-sumx/sumw1
		sig1<-sumxx/n
		dt<-max(c(max(abs(mu1-mu0)), max(abs(sig1-sig0))));
		mu0<-mu1;
		sig0<-sig1;
		n_it<-n_it+1;
	} ## end while
	if (n_it == max_it) warnings('Maximum iteration for EM algorithm is exceeded!')
	temp.cov<-sig1
	rownames(temp.cov)<-colnames(temp.cov)<-colnames(x)
	m<-matrix(mu1, p, 1)
	if (moment){
		p1<-p+1
		temp.cov<-matrix(NA, p+1,p+1)
		temp.cov[1,1]<-1
		temp.cov[1,2:p1]<-mu1
		temp.cov[2:p1, 1]<-mu1
		temp.cov[2:p1, 2:p1]<-sig1+m%x%t(m)
		rownames(temp.cov)<-colnames(temp.cov)<-c('(i)', colnames(x))
	}
	temp.cov
}

bmem.em<-function(x, ram, indirect, v, robust=FALSE, varphi=.1, st='i', moment=FALSE, max_it=500, ...){
	N<-dim(x)[1]
	if (missing(v)) v<-1:(ncol(x))
	temp.pattern<-bmem.pattern(x)
	if (robust){
		temp.cov<-bmem.em.rcov(temp.pattern, varphi, moment, max_it, st)
	}else{
		temp.cov<-bmem.em.cov(temp.pattern, moment, max_it)
	}
	temp.cov<-bmem.v(temp.cov,v,moment)
	bmem.sem(temp.cov, ram, N, indirect, moment, ...)
	}

bmem.list.jack<-function(x, ram, indirect, moment=FALSE, ...){
	n<-dim(x)[1]
	nseq<-1:n
	jack.est<-NULL
	jack.fit<-NULL
	for (i in 1:n){
		x.jack<-x[nseq[-i],]
		jack.temp<-try(bmem.list(x.jack, ram, indirect, moment, ...))
		if (class(jack.temp)!="try-error"){
			jack.est<-rbind(jack.est, jack.temp$est)
			jack.fit<-rbind(jack.fit, jack.temp$model.fit)
			}
		}
	list(jack.est=jack.est, jack.fit=jack.fit)
	}

bmem.pair.jack<-function(x, ram, indirect, moment=FALSE, ...){
	n<-dim(x)[1]
	nseq<-1:n
	jack.est<-NULL
	jack.fit<-NULL
	for (i in 1:n){
		x.jack<-x[nseq[-i],]
		jack.temp<-try(bmem.pair(x.jack, ram, indirect, moment, ...))
		if (class(jack.temp)!="try-error"){
			jack.est<-rbind(jack.est, jack.temp$est)
			jack.fit<-rbind(jack.fit, jack.temp$model.fit)
			}
		}
	list(jack.est=jack.est, jack.fit=jack.fit)
	}


bmem.mi.jack<-function(x, ram, indirect, v, m=10, moment=FALSE, ...){
	n<-dim(x)[1]
	nseq<-1:n
	jack.est<-NULL
	jack.fit<-NULL
	for (i in 1:n){
		x.jack<-x[nseq[-i],]
		jack.temp<-try(bmem.mi(x.jack, ram, indirect, v, m, moment, ...))
		if (class(jack.temp)!="try-error"){
			jack.est<-rbind(jack.est, jack.temp$est)
			jack.fit<-rbind(jack.fit, jack.temp$model.fit)
			}
		}
	list(jack.est=jack.est, jack.fit=jack.fit)
	}


bmem.em.jack<-function(x, ram, indirect, v, robust=FALSE, varphi=.1, st='i', moment=FALSE, max_it=500, ...){
	n<-dim(x)[1]
	nseq<-1:n
	jack.est<-NULL
	jack.fit<-NULL
	for (i in 1:n){
		x.jack<-x[nseq[-i],]
		jack.temp<-try(bmem.em(x.jack, ram, indirect, v, robust, varphi, st, moment, max_it, ...))
		if (class(jack.temp)!="try-error"){
			jack.est<-rbind(jack.est, jack.temp$est)
			jack.fit<-rbind(jack.fit, jack.temp$model.fit)
			}
		}
	list(jack.est=jack.est, jack.fit=jack.fit)
	}

bmem.list.boot<-function(x, ram, indirect, boot=1000, moment=FALSE, ...){	
	model0<-bmem.list(x, ram, indirect, moment, ...)
	par0<-model0$est
	fit0<-model0$model.fit
	n<-dim(x)[1]
	boot.est<-NULL
	boot.fit<-NULL
	for (i in 1:boot){
		x.boot<-x[sample(n,n, replace=TRUE),]
		modelb<-try(bmem.list(x.boot, ram, indirect, moment, ...))
		if (class(modelb)!="try-error"){
			boot.est<-rbind(boot.est, modelb$est)
			boot.fit<-rbind(boot.fit, modelb$model.fit)
			}
		}
	colnames(boot.fit)<-c('chisq', 'GFI','AGFI', 'RMSEA','NFI','NNFI','CFI','BIC','SRMR')	
	list(par.boot=boot.est, par0=par0, boot.fit=boot.fit, fit0=fit0)
	}

bmem.pair.boot<-function(x, ram, indirect, boot=1000, moment=FALSE, ...){
	model0<-bmem.pair(x, ram, indirect, moment, ...)
	par0<-model0$est
	fit0<-model0$model.fit
	n<-dim(x)[1]
	boot.est<-NULL
	boot.fit<-NULL
	for (i in 1:boot){
		x.boot<-x[sample(n, n, replace=TRUE),]
		modelb<-try(bmem.pair(x.boot, ram, indirect, moment,  ...))
		if (class(modelb)!="try-error"){
			boot.est<-rbind(boot.est, modelb$est)
			boot.fit<-rbind(boot.fit, modelb$model.fit)
			}
		}
	colnames(boot.fit)<-c('chisq', 'GFI','AGFI', 'RMSEA','NFI','NNFI','CFI','BIC','SRMR')	
	list(par.boot=boot.est, par0=par0, boot.fit=boot.fit, fit0=fit0)
	}

bmem.mi.boot<-function(x, ram, indirect, v, m=10, boot=1000, moment=FALSE, ...){
	model0<-bmem.mi(x, ram, indirect, v, m, moment, ...)
	par0<-model0$est
	fit0<-model0$model.fit
	n<-dim(x)[1]
	boot.est<-NULL
	boot.fit<-NULL
	for (i in 1:boot){
		x.boot<-x[sample(n, n, replace=TRUE),]
		modelb<-try(bmem.mi(x.boot, ram, indirect, v, m, moment, ...))
		if (class(modelb)!="try-error"){
			boot.est<-rbind(boot.est, modelb$est)
			boot.fit<-rbind(boot.fit, modelb$model.fit)
			}
		}
	colnames(boot.fit)<-c('chisq', 'GFI','AGFI', 'RMSEA','NFI','NNFI','CFI','BIC','SRMR')		
	list(par.boot=boot.est, par0=par0, boot.fit=boot.fit, fit0=fit0)
	}


bmem.em.boot<-function(x, ram, indirect, v, robust=FALSE, varphi=.1, st='i', boot=1000, moment=FALSE, max_it=500, ...){
	model0<-bmem.em(x, ram, indirect, v, robust, varphi, st, moment, max_it, ...)
	par0<-model0$est
	fit0<-model0$model.fit
	n<-dim(x)[1]
	boot.est<-NULL
	boot.fit<-NULL
	for (i in 1:boot){
		x.boot<-x[sample(n, n, replace=TRUE),]
		modelb<-try(bmem.em(x.boot, ram, indirect, v, robust, varphi, st, moment, max_it, ...))
		if (class(modelb)!="try-error"){
			boot.est<-rbind(boot.est, modelb$est)
			boot.fit<-rbind(boot.fit, modelb$model.fit)
			}
		}
	colnames(boot.fit)<-c('chisq', 'GFI','AGFI', 'RMSEA','NFI','NNFI','CFI','BIC','SRMR')	
	list(par.boot=boot.est, par0=par0, boot.fit=boot.fit, fit0=fit0)
	}
	
## bootstrap confidence intervals
bmem.ci.p<-function(par.boot, par0, cl=.95){
	alpha<-(1-cl)/2
	alpha<-c(alpha, 1-alpha)
	alpha<-sort(alpha)
	ci<-apply(par.boot, 2, quantile, prob=alpha, na.rm=TRUE)
	se.boot<-apply(par.boot, 2, sd, na.rm=TRUE)
	estimate<-par0
	cbind(estimate, se.boot, t(ci))
	}

bmem.ci.norm<-function(par.boot, par0, cl=.95){
	alpha<-(1-cl)/2
	alpha<-c(alpha, 1-alpha)
	alpha<-sort(alpha)
	se.boot<-apply(par.boot, 2, sd, na.rm=TRUE)
	estimate<-par0
	ci<-estimate+qnorm(alpha)%*%t(se.boot)
	dig <- max(2L, getOption("digits"))
	np<-length(alpha)
	qs <- paste(if (np < 100) 
            formatC(100 * alpha, format = "fg", width = 1, digits = dig)
        else format(100 * alpha, trim = TRUE, digits = dig), 
            "%", sep = "")
	ci.res<-cbind(estimate, se.boot, t(ci))
	colnames(ci.res)<-c('estimate','se.boot', qs)
	ci.res
	}

bmem.ci.bc1<-function(x, b, cl=.95){
	n<-length(x)
	z0<-qnorm(sum(x<b, na.rm=TRUE)/n)
	alpha<-(1-cl)/2
	alpha<-c(alpha, 1-alpha)
	alpha<-sort(alpha)
	alpha1<-alpha
	alpha<-pnorm(2*z0+qnorm(alpha))
	dig <- max(2L, getOption("digits"))
	np<-length(alpha)
	qs<-quantile(x, alpha, na.rm=TRUE)
	names(qs) <- paste(if (np < 100) 
            formatC(100 * alpha1, format = "fg", width = 1, digits = dig)
        else format(100 * alpha1, trim = TRUE, digits = dig), 
            "%", sep = "")
	qs
	}

bmem.ci.bc<-function(par.boot, par0, cl=.95){
	se.boot<-apply(par.boot, 2, sd, na.rm=TRUE)
	estimate<-par0
	p<-ncol(par.boot)
	ci<-NULL
	for (i in 1:p){
		ci<-rbind(ci, bmem.ci.bc1(par.boot[,i], par0[i], cl))
		}
	cbind(estimate, se.boot, ci)
	}

bmem.ci.bca1<-function(x, b, jack, cl=.95){
	n<-length(x)
	z0<-qnorm(sum(x<b, na.rm=TRUE)/n)
	alpha<-(1-cl)/2
	alpha<-c(alpha, 1-alpha)
	alpha<-sort(alpha)
	alpha1<-alpha
	mjack<-mean(jack, na.rm=TRUE)
	a<-sum((mjack-jack)^3, na.rm=TRUE)/(6*(sum((mjack-jack)^2, na.rm=TRUE))^1.5)
	alpha<-pnorm(z0+(z0+qnorm(alpha))/(1-a*(z0+qnorm(alpha))))
	dig <- max(2L, getOption("digits"))
	np<-length(alpha)
	if (alpha[1]=='NaN'){
		qs<-c(NA,NA)
	}else{
		qs<-quantile(x, alpha, na.rm=TRUE)
	}
	names(qs) <- paste(if (np < 100) 
            formatC(100 * alpha1, format = "fg", width = 1, digits = dig)
        else format(100 * alpha1, trim = TRUE, digits = dig), 
            "%", sep = "")
	qs
	}

bmem.ci.bca<-function(par.boot, par0, jack, cl=.95){
	se.boot<-apply(par.boot, 2, sd)
	estimate<-par0
	p<-ncol(par.boot)
	ci<-NULL
	for (i in 1:p){
		ci<-rbind(ci, bmem.ci.bca1(par.boot[,i], par0[i], jack[,i], cl))
		}
	cbind(estimate, se.boot, ci)
	}


### A main function for analysis
bmem<-function(x, ram, indirect, v, method='tsml', ci='bc', cl=.95, boot=1000, m=10, varphi=.1, st='i', robust=FALSE, max_it=500, moment=FALSE, ...){
	## method: list, pair, mi, me
	## ci: norm, perc, bc, bca
	N<-nrow(x)
	P<-ncol(x)
	if (missing(v)) v<-1:P
	
	## for listwise deletion method
	if (method=='list'){
		x<-x[,v]
		boot.est<-bmem.list.boot(x, ram, indirect, boot, moment, ...)
		if (ci=='norm'){
			ci.est<-bmem.ci.norm(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.norm(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='perc'){ 
			ci.est<-bmem.ci.p(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.p(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='bc'){ 
			ci.est<-bmem.ci.bc(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.bc(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='bca'){
			jack.est<-bmem.list.jack(x,ram,indirect,moment,...)
			ci.est<-bmem.ci.bca(boot.est$par.boot, boot.est$par0, jack.est$jack.est, cl)
			ci.fit<-bmem.ci.bca(boot.est$boot.fit, boot.est$fit0, jack.est$jack.fit, cl)
		}
	}
	
	## for pairwise deletion method
	if (method=='pair'){
		x<-x[,v]
		boot.est<-bmem.pair.boot(x, ram, indirect, boot, moment, ...)
		if (ci=='norm'){
			ci.est<-bmem.ci.norm(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.norm(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='perc'){ 
			ci.est<-bmem.ci.p(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.p(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='bc'){ 
			ci.est<-bmem.ci.bc(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.bc(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='bca'){
			jack.est<-bmem.pair.jack(x,ram,indirect,moment,...)
			ci.est<-bmem.ci.bca(boot.est$par.boot, boot.est$par0, jack.est$jack.est, cl)
			ci.fit<-bmem.ci.bca(boot.est$boot.fit, boot.est$fit0, jack.est$jack.fit, cl)
		}
	}
	
	## for multiple imputation method
	if (method=='mi'){
		boot.est<-bmem.mi.boot(x,ram,indirect,v,m,boot,moment,...)
		if (ci=='norm'){
			ci.est<-bmem.ci.norm(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.norm(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='perc'){ 
			ci.est<-bmem.ci.p(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.p(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='bc'){ 
			ci.est<-bmem.ci.bc(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.bc(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='bca'){
			jack.est<-bmem.mi.jack(x,ram,indirect,v,m,moment,...)
			ci.est<-bmem.ci.bca(boot.est$par.boot, boot.est$par0, jack.est$jack.est, cl)
			ci.fit<-bmem.ci.bca(boot.est$boot.fit, boot.est$fit0, jack.est$jack.fit, cl)
		}
	}
	
	## for EM method
	if (method=='tsml'){
		boot.est<-bmem.em.boot(x,ram, indirect, v, robust, varphi, st, boot, moment, max_it, ...)
		if (ci=='norm'){
			ci.est<-bmem.ci.norm(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.norm(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='perc'){ 
			ci.est<-bmem.ci.p(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.p(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='bc'){ 
			ci.est<-bmem.ci.bc(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.bc(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='bca'){
			jack.est<-bmem.em.jack(x,ram, indirect, v, robust, varphi, st, moment, max_it, ...)
			ci.est<-bmem.ci.bca(boot.est$par.boot, boot.est$par0, jack.est$jack.est, cl)
			ci.fit<-bmem.ci.bca(boot.est$boot.fit, boot.est$fit0, jack.est$jack.fit, cl)
		}
	}
	cat('The bootstrap confidence intervals for parameter estimates\n')
	print(ci.est)
	
	cat('\nThe bootstrap confidence intervals for model fit indices\n')
	rownames(ci.fit)<-c('chisq', 'GFI','AGFI', 'RMSEA','NFI','NNFI','CFI','BIC','SRMR')
	print(ci.fit)
	cat('\nThe literature has suggested the use of Bollen-Stine bootstrap for model fit. To do so, use the function bmem.bs().\n')
	if (ci=='bca') {
		bmemobject<-list(ci=ci.est, ci.fit=ci.fit, boot.est=boot.est, jack.est=jack.est)
	}else{
		bmemobject<-list(ci=ci.est, ci.fit=ci.fit, boot.est=boot.est)
	}
	class(bmemobject)<-'bmem'
	invisible(bmemobject)
}

summary.bmem<-function(object, ci='bc', cl=.95, ...){
	boot.est<-object$boot.est
		if (ci=='norm'){
			ci.est<-bmem.ci.norm(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.norm(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='perc'){ 
			ci.est<-bmem.ci.p(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.p(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='bc'){ 
			ci.est<-bmem.ci.bc(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.bc(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='bca'){
			jack.est<-object$jack.est
			ci.est<-bmem.ci.bca(boot.est$par.boot, boot.est$par0, jack.est$jack.est, cl)
			ci.fit<-bmem.ci.bca(boot.est$boot.fit, boot.est$fit0, jack.est$jack.fit, cl)
		}
	
	cat('The bootstrap confidence intervals for parameter estimates\n')
	print(ci.est)
	
	cat('\nThe bootstrap confidence intervals for model fit indices\n')
	rownames(ci.fit)<-c('chisq', 'GFI','AGFI', 'RMSEA','NFI','NNFI','CFI','BIC','SRMR')
	print(ci.fit)
	allci<-list(ci.est=ci.est, ci.fit=ci.fit)
	class(allci)<-'summary.bmem'
	invisible(allci)
}

bmem.raw2cov<-function(x){
	p<-nrow(x)
	m<-as.matrix(x[1, 2:p])
	Sraw<-x[2:p, 2:p]
	S<-Sraw-m%*%t(m)
	list(S=S, m=m)
	}

bmem.bs<-function(x, ram, indirect, v, ci='bc', cl=.95, boot=1000, max_it=500, ...){
	## Estimate the saturated mean and covariance matrix
	moment<-TRUE
	xmiss<-bmem.pattern(x)
	s.cov<-bmem.em.cov(xmiss, moment=TRUE, max_it)
	
	## Estimate the model indicated covariance matrix
	N<-nrow(x)
	temp.sem<-sem(ram, s.cov, N, raw=TRUE, fixed.x='(i)', ...)
	m.cov<-temp.sem$C
	
	## Take out the variables in the model (excluding the auxiliary variables)
	obsvar<-rownames(m.cov)
	s.cov.m<-s.cov[obsvar, obsvar]
	
	## mean and covariance matrix
	s.mcov<-bmem.raw2cov(s.cov.m)
	m.mcov<-bmem.raw2cov(m.cov)
	
	## select data
	x.m<-xmiss$x[,obsvar[2:length(obsvar)]]
	
	## generate new data based on the model
	x.bs<-x.m
	x.m.miss<-bmem.pattern(x.m)
	
	npat<-nrow(x.m.miss$misinfo)
	snj<-0
	for (j in 1:npat){
		nj<-x.m.miss$misinfo[j,1]
		pj<-x.m.miss$misinfo[j,2]
		oj<-x.m.miss$misinfo[j, 3:(2+pj)]
		mu_m<-m.mcov$m[oj,1]
		sig_m<-m.mcov$S[oj,oj]
		mu_s<-s.mcov$m[oj,1]
		sig_s<-s.mcov$S[oj,oj]
		for (i in (snj+1):(snj+nj)){
			x.bs[i, oj]<-(x.m.miss$x[i,oj] - t(mu_s)%*%solve(sig_s)%*%sig_m+t(mu_m))
			}
		snj<-nj	
		}
	
	## combine model data with auxiliary variables
	obsvar<-obsvar[-1]
	x[, obsvar]<-x.bs
	
	## now for bootstrap based on the newly generated sample
	N<-nrow(x)
	P<-ncol(x)
	if (missing(v)) v<-1:P
	boot.est<-bmem.em.boot(x,ram,indirect,v,boot,moment,max_it, ...)
		if (ci=='norm'){
			ci.est<-bmem.ci.norm(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.norm(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='perc'){ 
			ci.est<-bmem.ci.p(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.p(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='bc'){ 
			ci.est<-bmem.ci.bc(boot.est$par.boot, boot.est$par0, cl)
			ci.fit<-bmem.ci.bc(boot.est$boot.fit, boot.est$fit0, cl)
		}
		if (ci=='bca'){
			jack.est<-bmem.em.jack(x,ram,indirect,v,moment,max_it, ...)
			ci.est<-bmem.ci.bca(boot.est$par.boot, boot.est$par0, jack.est$jack.est, cl)
			ci.fit<-bmem.ci.bca(boot.est$boot.fit, boot.est$fit0, jack.est$jack.fit, cl)
		}
	cat('The bootstrap confidence intervals for parameter estimates\n')
	print(ci.est)
	
	cat('\nThe bootstrap confidence intervals for model fit indices\n')
	rownames(ci.fit)<-rownames(ci.fit)<-c('chisq', 'GFI','AGFI', 'RMSEA','NFI','NNFI','CFI','BIC','SRMR')
	print(ci.fit)
	cat('\nThe literature has suggested the use of raw data bootstrap for standard errors. To do so, use the function bmem().\n')

	if (ci=='bca') {
		bmemobject<-list(ci=ci.est, ci.fit=ci.fit, boot.est=boot.est, jack.est=jack.est)
	}else{
		bmemobject<-list(ci=ci.est, ci.fit=ci.fit, boot.est=boot.est)
	}
	class(bmemobject)<-'bmem'
	invisible(bmemobject)
}

bmem.plot<-function(x, par,...){
	## x: output from bmem function
	## par: a parameter or a fit indice to use
	parnames<-colnames(x$boot.est$par.boot)
	fitnames<-colnames(x$boot.est$boot.fit)
	
	if (par %in% parnames){
		hist(x$boot.est$par.boot[,par], xlab=par, ylab='Density', prob=TRUE, main='')
		lines(density(x$boot.est$par.boot[,par]))
	}else{
		hist(x$boot.est$boot.fit[,par], xlab=par, ylab='Density', prob=TRUE,main='')
		lines(density(x$boot.est$boot.fit[,par]))
	}	
}

plot.bmem<-function(x, par, ...){
	## x: output from bmem function
	## par: a parameter or a fit indice to use
	parnames<-colnames(x$boot.est$par.boot)
	fitnames<-colnames(x$boot.est$boot.fit)
	
	if (par %in% parnames){
		hist(x$boot.est$par.boot[,par], xlab=par, ylab='Density', prob=TRUE, main='')
		lines(density(x$boot.est$par.boot[,par]))
	}else{
		hist(x$boot.est$boot.fit[,par], xlab=par, ylab='Density', prob=TRUE,main='')
		lines(density(x$boot.est$boot.fit[,par]))
	}	
}

bmem.cov <- function(ram,obs.variables,moment=FALSE, debug=FALSE){
	if (moment){
		fixed.x<-'(i)'
	}else{
		fixed.x=NULL
	}
    parse.path <- function(path) {                                           
        path.1 <- gsub('-', '', gsub(' ','', path))
        direction <- if (regexpr('<>', path.1) > 0) 2 
            else if (regexpr('<', path.1) > 0) -1
            else if (regexpr('>', path.1) > 0) 1
            else stop(paste('ill-formed path:', path))
        path.1 <- strsplit(path.1, '[<>]')[[1]]
        list(first=path.1[1], second=path.1[length(path.1)], direction=direction)
        }
    if ((!is.matrix(ram)) | ncol(ram) != 3) stop ('ram argument must be a 3-column matrix')
    startvalues <- as.numeric(ram[,3])
    par.names <- ram[,2]
    n.paths <- length(par.names)
    heads <- from <- to <- rep(0, n.paths)
    for (p in 1:n.paths){
        path <- parse.path(ram[p,1])
        heads[p] <- abs(path$direction)
        to[p] <- path$second
        from[p] <- path$first
        if (path$direction == -1) {
            to[p] <- path$first
            from[p] <- path$second
            }
        }
    ram <- matrix(0, p, 5)
    all.vars <- unique(c(to, from))
    latent.vars <- setdiff(all.vars, obs.variables)
    not.used <- setdiff(obs.variables, all.vars)
    
    vars <- c(obs.variables, latent.vars)
    pars <- na.omit(unique(par.names))
    ram[,1] <- heads
    ram[,2] <- apply(outer(vars, to, '=='), 2, which)
    ram[,3] <- apply(outer(vars, from, '=='), 2, which)   
    par.nos <- apply(outer(pars, par.names, '=='), 2, which)
    if (length(par.nos) > 0)
        ram[,4] <- unlist(lapply(par.nos, function(x) if (length(x) == 0) 0 else x))
    ram[,5]<- startvalues
    colnames(ram) <- c('heads', 'to', 'from', 'parameter', 'start')
    if (!is.null(fixed.x)) fixed.x <- apply(outer(vars, fixed.x, '=='), 2, which)
    n <- length(obs.variables)
    m <- length(all.vars)
    t <- length(pars)
    if (debug) {
        cat('\n observed variables:\n') 
        print(paste(paste(1:n,':', sep=''), obs.variables, sep=''))
        cat('\n')
        if (m > n){ 
            cat('\n latent variables:\n')
            print(paste(paste((n+1):m,':', sep=''), latent.vars, sep=''))
            cat('\n')
            }
        cat('\n parameters:\n') 
        print(paste(paste(1:t,':', sep=''), pars, sep=''))
        cat('\n\n RAM:\n')
        print(ram)
        }
    n<-length(obs.variables)
    observed <- 1:n
    n.fix <- length(fixed.x)
    
    m <- max(ram[,c(2,3)])
    missing.variances <- setdiff(1:m, ram[,2][ram[,2] == ram[,3]])
    if (length(missing.variances) > 0) warning(paste( "\nThe model is almost surely misspecified; check also for missing covariances.\n"))
    t <- max(ram[,4])
    df <- n*(n + 1)/2 - t - n.fix*(n.fix + 1)/2
    if (df < 0) stop(paste("The model has negative degrees of freedom =", df))
    J <- matrix(0, n, m)
    correct <- matrix(2, m, m)
    diag(correct) <- 1
    J[cbind(1:n, observed)]<-1
    par.posn <-  sapply(1:t, function(i) which(ram[,4] == i)[1])
    colnames(ram)<-c("heads", "to", "from", "parameter", "start value")
    rownames(ram)<-rep("",nrow(ram))
    #if (length(param.names) > 0) rownames(ram)[par.posn]<-param.names
    fixed <- ram[,4] == 0
    sel.free <- ram[,4]
    sel.free[fixed] <- 1
    one.head <- ram[,1] == 1
    one.free <- which( (!fixed) & one.head )
    two.free <- which( (!fixed) & (!one.head) )
    arrows.1 <- ram[one.head, c(2,3), drop=FALSE]
    arrows.2 <- ram[!one.head, c(2,3), drop=FALSE]
    arrows.2t <- ram[!one.head, c(3,2), drop=FALSE]
    arrows.1.free <- ram[one.free,c(2,3), drop=FALSE]
    arrows.2.free <- ram[two.free,c(2,3), drop=FALSE]
    sel.free.1 <- sel.free[one.free]
    sel.free.2 <- sel.free[two.free]
    unique.free.1 <- unique(sel.free.1)
    unique.free.2 <- unique(sel.free.2)
    
    A <- P <- matrix(0, m, m)
    val <- ifelse (fixed, ram[,5], 0)
    A[arrows.1] <- val[one.head]
    P[arrows.2t] <- P[arrows.2] <- val[!one.head]
    I.Ainv <- solve(diag(m) - A)
    C <- J %*% I.Ainv %*% P %*% t(I.Ainv) %*% t(J)
    rownames(C)<-colnames(C)<-obs.variables
    C
    }

  