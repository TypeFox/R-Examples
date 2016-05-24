#-----------------------------------------------------------------------#
# Package: SAM      	   		   									    #
# Method: Sparse Additive Modelling using Quadratic Loss				#
# Authors: Tuo Zhao, Xingguo Li, Han Liu, and Kathryn Roeder            #
# Date: Mar 25th 2013                                                   #
# Version: 1.0.3                                                    	#
#-----------------------------------------------------------------------#

samQL = function(X, y, p=3, lambda = NULL, nlambda = NULL, lambda.min.ratio = 5e-3, thol=1e-5, max.ite = 1e5){
	
	gcinfo(FALSE)
	fit = list()
	fit$p = p
	
	fit = list()
	fit$p = p

	X = as.matrix(X)
	y = as.vector(y)

	n = nrow(X)
	d = ncol(X)
	m = d*p

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
	
	Z.mean = apply(Z,2,mean)
	Z.mean.rep = matrix(rep(Z.mean,n),nrow=n,byrow=T)
	Z = Z - Z.mean.rep
	y.mean = mean(y)
	y = y - y.mean	
	
	lambda_input = 1		
	if(is.null(lambda)){
		lambda_input = 0
		if(is.null(nlambda))
			nlambda = 30
		lambda = exp(seq(log(1),log(lambda.min.ratio),length = nlambda))
	} else nlambda = length(lambda)
	
	out = .C("grplasso",y = as.double(y), X = as.double(Z), lambda = as.double(lambda), nnlambda = as.integer(nlambda), nn = as.integer(n), dd = as.integer(d), pp = as.integer(p), ww = as.double(matrix(0,m,nlambda)), mmax_ite = as.integer(max.ite), tthol = as.double(thol),iinput = as.integer(lambda_input), df=as.integer(rep(0,nlambda)), sse=as.double(rep(0,nlambda)), func_norm = as.double(matrix(0,d,nlambda)), package="SAM")

	fit$lambda = out$lambda
	fit$w = matrix(out$w,ncol=nlambda)
	fit$df = out$df
	fit$sse = out$sse
	fit$func_norm = matrix(out$func_norm,ncol=nlambda)
	fit$intercept = rep(y.mean,nlambda) - t(Z.mean)%*%fit$w

	rm(out,X,y,Z,X.min.rep,X.ran.rep,Z.mean.rep)

	class(fit) = "samQL"
	return(fit)		
}

print.samQL = function(x,...){
	cat("Path length:",length(x$df),"\n")
	cat("d.f.:",x$df[1],"--->",x$df[length(x$df)],"\n")
}

plot.samQL = function(x,...){
	par = par(omi = c(0.0, 0.0, 0, 0), mai = c(1, 1, 0.1, 0.1))
	matplot(x$lambda[length(x$lambda):1],t(x$func_norm),type="l",xlab="Regularization Parameters",ylab = "Funcional Norms",cex.lab=2,log="x",lwd=2)
}

predict.samQL = function(object, newdata,...){
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
	
	out$values = cbind(Zt,rep(1,nt))%*%rbind(object$w,object$intercept)
	
	rm(Zt,newdata)
	
	return(out)
}

