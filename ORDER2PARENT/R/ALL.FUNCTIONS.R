bgmm<-function(dat,orderinfo,degree = 3,support = NULL,weight.type = 1){
    # Implementing the B-spline GMM estimator.
    # dat: a list of observations on various order statistics.
    # orderinfo: a J*2 matrix containing the rank as well as the sample size of all J order statistics.
    # degree: the degree of B-spline. The default value is 3, i.e. cubic B-spline.
    # support: the support of the parent distribution. The default is NULL. If null, the minimum and maximum of the sample will be used as the support.
    # weight.type: which type of weight matrix the user is going to use. The first type is based on sample size, and the second one is based on initial estimation MSE.
    library(splines)
    library(Matrix)
    if(nrow(orderinfo)!=length(dat)){
	stop('The information about order statistics is missing')
    }
    J<-nrow(orderinfo)
    
    # Compute the ECDF functions for each order statistics.
    ecdf.order<-list()	# a list containing the ECDF function for each order statistic.
    ecdf.parent<-list()	# a list containing the parent ECDF function for each order statistic.
    ecdf.parent.vector<-c()
    x<-c()	# the whole data.
    n.order<-c()
    for(j in 1:J){
	n.order[j]<-length(dat[[j]])
	ecdf.order[[j]]<-ecdf(dat[[j]])
	ecdf.parent[[j]]<-parentcdf(ecdf.order[[j]](dat[[j]]),k=orderinfo[j,1],m=orderinfo[j,2])
	ecdf.parent.vector<-c(ecdf.parent.vector,ecdf.parent[[j]])
	x<-c(x,dat[[j]])
    }
    n<-length(x)
    
    # Each knots selection, implement the following codes.
    SplineGMM<-function(EF, basis.x, type = 1, weight.matrix = NULL){
	# the default weight.matrix is to use the sample size based matrix.
	# type control the type of weight matrix in use.
	moment.conditions<-function(beta.hat){
	    # beta.hat is the estimator.
	    F0<-basis.x%*%beta.hat	# a n*1 vector.
	    end.i<-c()
	    for(j in 1:J){
		end.i[j]<-sum(n.order[1:j])
	    }
	    start.i<-end.i-n.order+1
	    n.para<-ncol(basis.x)
	    sample.moments<-matrix(nr=n.para,nc=J)
	    for(j in 1:J){
		dev<-EF[[j]]-F0[start.i[j]:end.i[j]]
		if(type==1){
		    sample.moments[,j]<-(t(basis.x[start.i[j]:end.i[j],])%*%dev)	# (q+p)*1 vector.
		}else{
		    sample.moments[,j]<-(t(basis.x[start.i[j]:end.i[j],])%*%dev)/n.order[j]	# (q+p)*1 vector.
		}
	    }
	    return(as.vector(sample.moments))
	}
	loss.gmm<-function(betahat,type,weight.matrix=NULL){
	    moments.betahat<-moment.conditions(betahat)
	    if(type==1){
		q<-t(moments.betahat)%*%moments.betahat
	    }else{
		q<-t(moments.betahat)%*%weight.matrix%*%moments.betahat
	    }
	    return(as.double(q))
	}
	n.para<-ncol(basis.x)
	constraint<-cbind(diag(-1,n.para-1),rep(0,n.para-1),deparse.level=0)+cbind(rep(0,n.para-1),diag(1,n.para-1),deparse.level=0)
	betahat<-constrOptim(seq(from=.1,to=1,length.out=n.para),f=loss.gmm,ui=constraint,ci=rep(0,n.para-1),method='Nelder-Mead',type=type,weight.matrix=weight.matrix)$par
	return(betahat)
    }
    
    weight.error<-function(EF,basis.x,betahat){
	# compute the pseudo-weight matrix.
	F0<-basis.x%*%betahat	# the estimates of parent cdf based on the first-stage estimates.
	end.i<-c()
	for(j in 1:J){
	    end.i[j]<-sum(n.order[1:j])
	}
	start.i<-end.i-n.order+1
	weight.m<-list()
	width.basis<-ncol(basis.x)
	for(j in 1:J){
	    w.m<-diag(mean((EF[[j]]-F0[start.i[j]:end.i[j]])^2),width.basis)
	    weight.m[[j]]<-w.m
	}
	weight.m<-solve(bdiag(weight.m))
	return(weight.m)
    }
    
    est.stage1<-blr(dat,orderinfo,degree,support,constraint=TRUE)
    n.knots<-est.stage1$n.knots
    betahat1<-est.stage1$betahat
    knots.bs<-c(1:n.knots)/(n.knots+1)
    knots.bs<-sort(x)[as.integer(n*knots.bs)]
    if(is.null(support)){
	v.lower<-min(x)
	v.upper<-max(x)
    }else{
	v.lower<-min(support)
	v.upper<-max(support)
    }
    basis.x<-bs(x,knots=knots.bs,degree=degree,Boundary.knots=c(v.lower,v.upper))	# basis matrix.
    if(weight.type==1){
	beta.gmm<-SplineGMM(ecdf.parent,basis.x)
    }else{
	# use error based weight matrix.
	weightinuse<-weight.error(ecdf.parent,basis.x,betahat1)
	beta.gmm<-SplineGMM(ecdf.parent,basis.x,type=weight.type,weight.matrix=weightinuse)
    }
    return(list(betahat=beta.gmm,n.knots=n.knots))
}

blr<-function(dat,orderinfo,degree = 3,support = NULL,constraint = FALSE){
    # Implementing the B-spline linear regression estimator.
    # dat: a list of observations on various order statistics.
    # orderinfo: a J*2 matrix containing the rank as well as the sample size of all J order statistics.
    # degree: the degree of B-spline. The default value is 3, i.e. cubic B-spline.
    # support: the support of the parent distribution. The default is NULL. If null, the minimum and maximum of the sample will be used as the support.
    # contraint: a logical variable indicating whether implement the contrainted B-spline linear regression estimator.
    library(splines)
    if(nrow(orderinfo)!=length(dat)){
	stop('The information about order statistics is missing')
    }
    J<-nrow(orderinfo)
    
    # Compute the ECDF functions for each order statistics.
    ecdf.order<-list()	# a list containing the ECDF function for each order statistic.
    ecdf.parent<-c()	# a list containing the parent ECDF function for each order statistic.
    x<-c()	# the whole data.
    for(j in 1:J){
	ecdf.order[[j]]<-ecdf(dat[[j]])
	ecdf.parent<-c(ecdf.parent,parentcdf(ecdf.order[[j]](dat[[j]]),k=orderinfo[j,1],m=orderinfo[j,2]))
	x<-c(x,dat[[j]])
    }
    n<-length(x)
    q.tilta<-as.integer(n^(1/(2*degree+3)))
    n.knots<-c(max(1,.5*q.tilta):min(5*q.tilta,n/4-degree))	# number of inner knots candidates.
    n.candidates<-length(n.knots)
    
    # square loss function
    squareloss<-function(betahat,basis){
	dev<-ecdf.parent-basis%*%betahat
	loss<-as.double(t(dev)%*%dev)/n
	return(loss)
    }
    beta.list<-list()
    AIC.bs<-c()
    for(k in 1:n.candidates){
	knots.bs<-c(1:n.knots[k])/(n.knots[k]+1)
	knots.bs<-sort(x)[as.integer(n*knots.bs)]
	if(is.null(support)){
	    v.lower<-min(x)
	    v.upper<-max(x)
	}else{
	    v.lower<-min(support)
	    v.upper<-max(support)
	}
	basis.x<-bs(x,knots=knots.bs,degree=degree,Boundary.knots=c(v.lower,v.upper))	# basis matrix.
	if(constraint){
	    np<-n.knots[k]+degree	# number of unknown parameters
	    constraint.m<-cbind(diag(-1,np-1),rep(0,np-1),deparse.level=0)+cbind(rep(0,np-1),diag(1,np-1),deparse.level=0)
	    blm<-constrOptim(seq(from=.1,to=1,length.out=np),f=squareloss,ui=constraint.m,ci=rep(0,np-1),method='Nelder-Mead',basis=basis.x)
	    beta.ls<-blm$par
	    loss.ls<-blm$value
	    beta.list[[k]]<-beta.ls
	}else{
	    blm<-lm(ecdf.parent~basis.x-1)
	    beta.ls<-blm$coef
	    loss.ls<-mean((blm$residuals)^2)
	    beta.list[[k]]<-beta.ls
	}
	AIC.bs[k]<-log(loss.ls)+2*(n.knots[k]+degree)/n
    }
    beta.ls<-beta.list[[which.min(AIC.bs)]]
    n.knots<-n.knots[which.min(AIC.bs)]
    return(list(betahat=beta.ls,n.knots=n.knots))
}

parentcdf<-function(F.order,k,m){
    # Transform the cdf of the k-th order statistic in random samples of m to its parent cdf.
    # The input F.order can either be a scalar or a vector.
    combinations<-function(x,m){
	# calculate combinations number.
	if(m==0){
	    return(1)
	}else{
	    return(ncol(combn(x,m)))
	}
    }
    RHStransformer<-function(F,k,m){
	# F is the parent distribution.
	# This function computes the RHS of the equation.
	n.box<-(k-1)+c(0:(m-k))
	v.combn<-sapply(n.box,combinations,m=(k-1))
	v.F<-(1-F)^c(0:(m-k))
	RHS<-(F^k)*(t(v.combn)%*%v.F)
	return(RHS)
    }
    solver<-function(F.order,k,m){
	# Currently, I only allow m=k+1 to avoid too much computation burden.
	dev<-function(F,F.order){
	    F.order-RHStransformer(F,k,m)
	}
	F.est<-uniroot(dev,interval=c(0,1),F.order=F.order)$root
	return(F.est)
    }
    F0<-sapply(F.order,solver,k=k,m=m)
    return(F0)
}


parentest<-function(x0,beta.hat,n.knots,degree = 3,support = NULL){
    # for any x, this function will return our estimate of its cdf.
    n0<-length(x0)
    knots.bs<-c(1:n.knots)/(n.knots+1)
    knots.bs<-sort(x0)[as.integer(n0*knots.bs)]
    if(is.null(support)){
	v.lower<-min(x0)
	v.upper<-max(x0)
    }else{
	v.lower<-min(support)
	v.upper<-max(support)
    }
    basis.x0<-bs(x0,knots=knots.bs,degree=degree,Boundary.knots=c(v.lower,v.upper))	# basis matrix.
    F0hat<-basis.x0%*%beta.hat
    return(F0hat)
}