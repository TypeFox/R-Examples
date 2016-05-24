# Functions from package 'SKAT' v.0.82 (c) 2011, with changes

#	Function get parameters of optimal test
#
SKAT_Optimal_Param<-function(Z1,r.all){

	n<-dim(Z1)[1]
	p.m<-dim(Z1)[2]
	r.n<-length(r.all)

	z_mean<-rowMeans(Z1)
	Z_mean<-matrix(rep(z_mean,p.m),ncol=p.m,byrow=FALSE)
	cof1<-(t(z_mean) %*% Z1)[1,] / sum(z_mean^2)

	Z.item1<-Z_mean %*% diag(cof1)	
	Z.item2<-Z1 - Z.item1

	# W3.2 Term : mixture chisq
	W3.2.t<-t(Z.item2) %*% Z.item2
	lambda<-Get_Lambda(W3.2.t)
	
	# W3.3 Term : variance of remaining ...
	W3.3.item<-sum((t(Z.item1) %*% Z.item1) * (t(Z.item2) %*% Z.item2)) * 4
	
	# Mixture Parameters
	MuQ<-sum(lambda)
	VarQ<-sum(lambda^2) *2 + W3.3.item
	KerQ<-sum(lambda^4)/(sum(lambda^2))^2 * 12
	Df<-12/KerQ

	# W3.1 Term : tau1 * chisq_1
	tau<-rep(0,r.n)
	for(i in 1:r.n){
		r.corr<-r.all[i]
		#term1<-p.m*r.corr + cof1^2 * (1-r.corr)
		term1<-p.m^2*r.corr + sum(cof1^2) * (1-r.corr)
		tau[i]<-sum(term1) *  sum(z_mean^2)
	}

	out<-list(MuQ=MuQ,VarQ=VarQ,KerQ=KerQ,lambda=lambda,VarRemain=W3.3.item,Df=Df,tau=tau)
	return(out)
}

#
#	Function get SKAT statistics with given rho
#		Q.all is a matrix with n.q x n.r

SKAT_Optiaml_Each_Q<-function(param.m, Q.all, r.all, lambda.all, method=NULL){

	n.r<-length(r.all)
	c1<-rep(0,4)
	n.q<-dim(Q.all)[1]

	pval<-matrix(rep(0,n.r*n.q),ncol=n.r)
	pmin.q<-matrix(rep(0,n.r*n.q),ncol=n.r)
	param.mat<-NULL

	for(i in 1:n.r){
		Q<-Q.all[,i]
		r.corr<-r.all[i]
		lambda.temp<-lambda.all[[i]] 
		c1[1]<-sum(lambda.temp)
		c1[2]<-sum(lambda.temp^2)
		c1[3]<-sum(lambda.temp^3)
		c1[4]<-sum(lambda.temp^4)
		param.temp<-Get_Liu_Params_Mod(c1)

		muQ<-param.temp$muQ
		varQ<-param.temp$sigmaQ^2
		df<-param.temp$l

		# get pvalue
		Q.Norm<-(Q - muQ)/sqrt(varQ) * sqrt(2*df) + df
		pval[,i]<- pchisq(Q.Norm,  df = df, lower.tail=FALSE)
		# will be changed later
		
		if(!is.null(method)){
			if(method=="optimal.mod" || method=="optimal.adj" || method=="optimal.moment.adj" ){
				pval[,i]<-Get_PValue.Lambda(lambda.temp,Q)$p.value
			}
		}
		
		param.mat<-rbind(param.mat,c(muQ,varQ,df))
	}

	pmin<-apply(pval,1,min)
	for(i in 1:n.r){
	
		muQ<-param.mat[i,1]
		varQ<-param.mat[i,2]
		df<-param.mat[i,3]

		q.org<-qchisq(1-pmin,df=df)
		q.q<-(q.org - df)/sqrt(2*df) *sqrt(varQ) + muQ
		pmin.q[,i]<-q.q

	}
	
	out<-list(pmin=pmin,pval=pval,pmin.q=pmin.q)
	return(out)

}

SKAT_Optimal_Integrate_Func_Davies<-function(x,pmin.q,param.m,r.all,acc,lim){
	
	n.r<-length(r.all)
	n.x<-length(x)

	temp1<-param.m$tau %x% t(x)

	temp<-(pmin.q - temp1)/(1-r.all)
	temp.min<-apply(temp,2,min)

	re<-rep(0,length(x))
	for(i in 1:length(x)){
		#a1<<-temp.min[i]
		min1<-temp.min[i]
		if(min1 > sum(param.m$lambda) * 10^4){
			temp<-0
		} else {
			min1.temp<- min1 - param.m$MuQ			
			sd1<-sqrt(param.m$VarQ - param.m$VarRemain)/sqrt(param.m$VarQ)
			min1.st<-min1.temp *sd1 + param.m$MuQ
			
			dav.re<-davies(min1.st,param.m$lambda,acc = acc,lim = lim)
			temp<-dav.re$Qq
			if(dav.re$ifault != 0){
				stop("dav.re$ifault is not 0")
			}
		}
		if(temp > 1){
			temp=1
		}
		#lambda.record<<-param.m$lambda
		#print(c(min1,temp,dav.re$ifault,sum(param.m$lambda)))
		re[i]<-(1-temp) * dchisq(x[i],df=1)
	}
	return(re)

}

SKAT_Optimal_Integrate_Func_Kuonen<-function(x,pmin.q,param.m,r.all){
	
	n.r<-length(r.all)
	n.x<-length(x)

	temp1<-param.m$tau %x% t(x)

	temp<-(pmin.q - temp1)/(1-r.all)
	temp.min<-apply(temp,2,min)

	re<-rep(0,length(x))
	for(i in 1:length(x)){
		#a1<<-temp.min[i]
		min1<-temp.min[i]
		if(min1 > sum(param.m$lambda) * 10^4){
			temp<-0
		} else {
			min1.temp<- min1 - param.m$MuQ			
			sd1<-sqrt(param.m$VarQ - param.m$VarRemain)/sqrt(param.m$VarQ)
			min1.st<-min1.temp *sd1 + param.m$MuQ
			temp <- pchisqsum(min1.st, rep(1, length(param.m$lambda)), param.m$lambda, lower.tail = F, method = 'sad')
		}
		if(temp > 1){
			temp=1
		}
		#lambda.record<<-param.m$lambda
		#print(c(min1,temp,dav.re$ifault,sum(param.m$lambda)))
		re[i]<-(1-temp) * dchisq(x[i],df=1)
	}
	return(re)

}

# add pmin on 02-13-2013 
SKAT_Optimal_PValue_Davies_Kuonen<-function(pmin.q,param.m,r.all, pmin=NULL, method, acc, lim){

	#re<-try(integrate(SKAT_Optimal_Integrate_Func_Davies, lower=0, upper=30, subdivisions=500, pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-15), silent = TRUE)

	if (method == 'kuonen') { re<-try(integrate(SKAT_Optimal_Integrate_Func_Kuonen, lower=0, upper=40, subdivisions=1000, pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-25), silent = TRUE)
	} else { re<-try(integrate(SKAT_Optimal_Integrate_Func_Davies, lower=0, upper=40, subdivisions=1000, pmin.q=pmin.q,param.m=param.m,r.all=r.all,acc=acc,lim=lim,abs.tol = 10^-25), silent = TRUE) }

	if(class(re) == "try-error"){
		re<-SKAT_Optimal_PValue_Liu(pmin.q,param.m,r.all, pmin)
		return(re)
	} 

	pvalue<-1-re[[1]]
	if(!is.null(pmin)){
		if(pmin *length(r.all) < pvalue){
			pvalue = pmin *length(r.all)
		}
	}
	
	
	return(pvalue)

}


SKAT_Optimal_Integrate_Func_Liu<-function(x,pmin.q,param.m,r.all){
	
	#x<-1
	#print(length(x))
	#print(x)
	#X1<<-x
	#x<-X1

	n.r<-length(r.all)
	n.x<-length(x)

	temp1<-param.m$tau %x% t(x)

	temp<-(pmin.q - temp1)/(1-r.all)
	temp.min<-apply(temp,2,min)

	temp.q<-(temp.min - param.m$MuQ)/sqrt(param.m$VarQ)*sqrt(2*param.m$Df) + param.m$Df
	re<-pchisq(temp.q ,df=param.m$Df) * dchisq(x,df=1)
	
	return(re)

}

# add pmin on 02-13-2013 
SKAT_Optimal_PValue_Liu<-function(pmin.q,param.m,r.all, pmin=NULL){

	 re<-integrate(SKAT_Optimal_Integrate_Func_Liu, lower=0, upper=40, subdivisions=2000
	,pmin.q=pmin.q,param.m=param.m,r.all=r.all,abs.tol = 10^-25)
	
	pvalue<-1-re[[1]]
	
	if(!is.null(pmin)){
		if(pmin *length(r.all) < pvalue){
			pvalue = pmin *length(r.all)
		}
	}
	
	return(pvalue)

}

SKAT_Optimal_Get_Pvalue<-function(Q.all, Z1, r.all, method, acc, lim){

	n.r<-length(r.all)
	n.q<-dim(Q.all)[1]
	p.m<-dim(Z1)[2]

	lambda.all<-list()
	for(i in 1:n.r){
		r.corr<-r.all[i]
		R.M<-diag(rep(1-r.corr,p.m)) + matrix(rep(r.corr,p.m*p.m),ncol=p.m)
		if (r.corr < 1) { L<-chol(R.M,pivot=TRUE)
		} else {
		L <- R.M / sqrt(p.m)
		}
		Z2<- Z1 %*% t(L)
		K1<-t(Z2) %*% Z2

		lambda.all[[i]]<-Get_Lambda(K1)
		 
	}

	# Get Mixture param 
	param.m<-SKAT_Optimal_Param(Z1,r.all)
	Each_Info<-SKAT_Optiaml_Each_Q(param.m, Q.all, r.all, lambda.all, method=method)
	pmin.q<-Each_Info$pmin.q
	pmin<-Each_Info$pmin
	pval<-rep(0,n.q)

	if(method == "davies" || method=="kuonen" || method=="optimal" || method=="optimal.mod" || method=="optimal.adj"){

		for(i in 1:n.q){
			pval[i]<-SKAT_Optimal_PValue_Davies_Kuonen(pmin.q[i,],param.m,r.all, pmin[i], method, acc = acc, lim = lim)
			
		}

	} else if(method =="liu" || method =="liu.mod" || method=="optimal.moment" || method=="optimal.moment.adj" ){
		
		for(i in 1:n.q){
			pval[i]<-SKAT_Optimal_PValue_Liu(pmin.q[i,],param.m,r.all, pmin[i])
		}

	} else {
		stop("Invalid Method!")
	}
	
	# Check the pval 
	# Since SKAT-O is between burden and SKAT, SKAT-O p-value should be <= min(p-values) * 2
	# To correct conservatively, we use min(p-values) * 3
	
	multi<-3
	if(length(r.all) < 3){
		multi<-2
	}

	for(i in 1:n.q){
		pval.each<-Each_Info$pval[i,]
		IDX<-which(pval.each > 0)
		
		pval1<-min(pval.each) * multi
		if(pval[i] <= 0 || length(IDX) < length(r.all)){
			pval[i]<-pval1
		}
		
		# if pval==0, use nonzero min each.pval as p-value
		if(pval[i] == 0){
			if(length(IDX) > 0){
				pval[i] = min(pval.each[IDX])
			}
		}
	
	}
	
	return(list(p.value=pval,p.val.each=Each_Info$pval))

}
