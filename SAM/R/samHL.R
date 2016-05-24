#-----------------------------------------------------------------------#
# Package: SAM      	   		   									    #
# Method: Sparse Additive Modelling using Hinge Loss					#
# Authors: Tuo Zhao, Xingguo Li, Han Liu, and Kathryn Roeder            #
# Date: Mar 25th 2013                                                   #
# Version: 1.0.3                                                    	#
#-----------------------------------------------------------------------#

samHL = function(X, y, p=3, lambda = NULL, nlambda = NULL, lambda.min.ratio = 0.4, thol=1e-5, mu = 5e-2, max.ite = 1e5){
	
	gcinfo(FALSE)
	fit = list()
	fit$p = p
	
	fit = list()
	
	X = as.matrix(X)
	y = as.vector(y)
	
	
	n = nrow(X)
	d = ncol(X)
	m = d*p

	np = sum(y==1)
	nn = sum(y==-1)

	if((np+nn)!=n){
		cat("Please check the labels. (Must be coded in 1 and -1)")
		fit = "Please check the labels."
		return(fit)
	}
	
	if(np>nn) a0 = 1-nn/np*mu else a0 = np/nn*mu - 1
	
	fit$p = p

	X.min = apply(X,2,min)
	X.max = apply(X,2,max)
	X.ran = X.max - X.min
	X.min.rep = matrix(rep(X.min,n),nrow=n,byrow=T)
	X.ran.rep = matrix(rep(X.ran,n),nrow=n,byrow=T)
	X = (X-X.min.rep)/X.ran.rep

	fit$X.min = X.min
	fit$X.ran = X.ran	
	
	Z = matrix(0,n,m)
	fit$nkots = matrix(0,p-1,d)
	fit$Boundary.knots = matrix(0,2,d)
	for(j in 1:d){
		tmp = (j-1)*p + c(1:p)
		tmp0 = ns(X[,j],df=p)
		Z[,tmp] = tmp0
		fit$nkots[,j] = attr(tmp0,'knots')
		fit$Boundary.knots[,j] = attr(tmp0,'Boundary.knots')
	}
	
	Z = cbind(matrix(rep(y,m),n,m)*Z,y)	
	
	if(is.null(lambda)){
		u = cbind((rep(1,n) - a0*y)/mu,rep(0,n),rep(1,n))
		u = apply(u,1,median)
		
		if(is.null(nlambda)) nlambda = 20;
		
		lambda_max = max(sqrt(colSums(matrix(t(Z[,1:(p*d)])%*%u,p,d)^2)))
		lambda = exp(seq(log(1),log(lambda.min.ratio),length=nlambda))*lambda_max
	}
	else nlambda = length(lambda)
	
	L0 = norm(Z,'f')^2/mu
	
	out = .C("grpSVM", Z = as.double(Z), lambda = as.double(lambda), nnlambda = as.integer(nlambda), LL0 = as.double(L0), nn = as.integer(n), dd = as.integer(d), pp = as.integer(p),aa0 = as.double(a0), xx = as.double(matrix(0,m+1,nlambda)), mmu = as.double(mu), mmax_ite = as.integer(max.ite), tthol = as.double(thol),aalpha = as.double(0.5),df=as.double(rep(0,nlambda)),func_norm=as.double(matrix(0,d,nlambda)),package="SAM")

	fit$lambda = out$lambda
	fit$w = matrix(out$xx,ncol=nlambda)
	fit$df = out$df
	fit$func_norm = matrix(out$func_norm,ncol=nlambda)
	
	rm(out,X,y,Z,X.min.rep,X.ran.rep)

	class(fit) = "samHL"
	return(fit)		
}

print.samHL = function(x,...){
	cat("Path length:",length(x$df),"\n")
	cat("d.f.:",x$df[1],"--->",x$df[length(x$df)],"\n")
}

plot.samHL = function(x,...){
	par = par(omi = c(0.0, 0.0, 0, 0), mai = c(1, 1, 0.1, 0.1))
	matplot(x$lambda[length(x$lambda):1],t(x$func_norm),type="l",xlab="Regularization Parameters",ylab = "Funcional Norms",cex.lab=2,log="x",lwd=2)
}

predict.samHL = function(object, newdata, thol = 0, ...){
	gcinfo(FALSE)
	out = list()
	nt = nrow(newdata)
	d = ncol(newdata)
	X.min.rep = matrix(rep(object$X.min,nt),nrow=nt,byrow=T)
	X.ran.rep = matrix(rep(object$X.ran,nt),nrow=nt,byrow=T)
		
	newdata = (newdata-X.min.rep)/X.ran.rep
	newdata = pmax(newdata,0)
	newdata = pmin(newdata,1)
	
	m = object$p*d
	
	Zt = matrix(0,nt,m)

	for(j in 1:d){
		tmp = (j-1)*object$p + c(1:object$p)
		Zt[,tmp] = ns(newdata[,j],df=object$p,knots=object$knots[,j],Boundary.knots=object$Boundary.knots[,j])
	}
	
	out$values = cbind(Zt,rep(1,nt))%*%object$w
	out$labels = sign(out$values-0)
	
	rm(Zt,newdata)
	
	return(out)
}

