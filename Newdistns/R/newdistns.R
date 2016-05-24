dgammag<-function(x, spec, a=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-f(x, ...)*(1-F(x, ...))**(-2)*dgamma(F(x, ...)/(1-F(x, ...)),shape=a)
        pdf[log==TRUE]<-fL(x, ...)-2*log(1-F(x, ...))+dgamma(F(x, ...)/(1-F(x, ...)),shape=a,log=TRUE)
        return(pdf)

}



pgammag<-function(x, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-pgamma(F(x, ...)/(1-F(x, ...)),shape=a)
        cdf[log.p==TRUE&lower.tail==TRUE]<-pgamma(F(x, ...)/(1-F(x, ...)),shape=a,log.p=TRUE)
        cdf[log.p==FALSE&lower.tail==FALSE]<-1-pgamma(F(x, ...)/(1-F(x, ...)),shape=a)
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(1-pgamma(F(x, ...)/(1-F(x, ...)),shape=a))
        return(cdf)

}



qgammag<-function(p, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq(qgamma(p[p>=0&p<=1],shape=a)/(1+qgamma(p[p>=0&p<=1],shape=a)), ...)
	return(qf)

}


rgammag<-function(n, spec, a=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qgammag(u, spec, a=a, ...)
	return(sf)
}




dbetaexpg<-function(x, spec, lambda=1, a=1, b=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-lambda*f(x, ...)*(1-F(x, ...))**(b-1)*dbeta(1-(1-F(x, ...))**lambda,shape1=a,shape2=b)
        pdf[log==TRUE]<-log(lambda)+fL(x, ...)+(b-1)*log(1-F(x, ...))+dbeta(1-(1-F(x, ...))**lambda,shape1=a,shape2=b,log=TRUE)
        return(pdf)

}


pbetaexpg<-function(x, spec, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-1-pbeta((1-F(x, ...))**lambda,shape1=lambda*(b-1)+1,shape2=a)
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(1-pbeta((1-F(x, ...))**lambda,shape1=lambda*(b-1)+1,shape2=a))
        cdf[log.p==FALSE&lower.tail==FALSE]<-pbeta((1-F(x, ...))**lambda,shape1=lambda*(b-1)+1,shape2=a)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pbeta((1-F(x, ...))**lambda,shape1=lambda*(b-1)+1,shape2=a,log.p=TRUE)
	return(cdf)

}


qbetaexpg<-function(p, spec, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq(1-(qbeta(1-p[p>=0&p<=1],shape1=lambda*(b-1)+1,shape2=a))**(1/lambda), ...)
	return(qf)

}


rbetaexpg<-function(n, spec, lambda=1, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qbetaexpg(u, spec, lambda=lambda, a=a, b=b, ...)
	return(sf)
}


dweibullg<-function(x, spec, beta=1, c=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-c*beta**(-c)*(f(x, ...)/(1-F(x, ...)))*(-log(1-F(x, ...)))**(c-1)*exp(-beta**(-c)*(-log(1-F(x, ...)))**c)
        pdf[log==TRUE]<-log(c)-c*log(beta)+fL(x, ...)-log(1-F(x, ...))+(c-1)*log(-log(1-F(x, ...)))-beta**(-c)*(-log(1-F(x, ...)))**c
        return(pdf)

}


pweibullg<-function(x, spec, beta=1, c=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-1-exp(-beta**(-c)*(-log(1-F(x, ...)))**c)
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(1-exp(-beta**(-c)*(-log(1-F(x, ...)))**c))
        cdf[log.p==FALSE&lower.tail==FALSE]<-exp(-beta**(-c)*(-log(1-F(x, ...)))**c)
        cdf[log.p==TRUE&lower.tail==FALSE]<--beta**(-c)*(-log(1-F(x, ...)))**c
        return(cdf)

}


qweibullg<-function(p, spec, beta=1, c=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq(1-exp(-beta*(-log(1-p[p>=0&p<=1]))**(1/c)), ...)
	return(qf)

}


rweibullg<-function(n, spec, beta=1, c=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qweibullg(u, spec, beta=beta, c=c, ...)
	return(sf)
}





dgepg<-function(x, spec, theta=1, eta=0.5, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-theta*(1-eta)*(1-exp(-theta))*f(x, ...)*exp(-theta+theta*F(x, ...))/(1-exp(-theta)-eta+eta*exp(-theta+theta*F(x, ...)))**2
        pdf[log==TRUE]<-log(theta)+log(1-eta)+log(1-exp(-theta))+fL(x, ...)-theta+theta*F(x, ...)-2*log(1-exp(-theta)-eta+eta*exp(-theta+theta*F(x, ...)))
	return(pdf)

}


pgepg<-function(x, spec, theta=1, eta=0.5, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-(exp(-theta+theta*F(x, ...))-exp(-theta))/(1-exp(-theta)-eta+eta*exp(-theta+theta*F(x, ...)))
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(exp(-theta+theta*F(x, ...))-exp(-theta))-log(1-exp(-theta)-eta+eta*exp(-theta+theta*F(x, ...)))
        cdf[log.p==FALSE&lower.tail==FALSE]<-(1-eta)*(1-exp(-theta+theta*F(x, ...)))/(1-exp(-theta)-eta+eta*exp(-theta+theta*F(x, ...)))
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(1-eta)+log(1-exp(-theta+theta*F(x, ...)))-log(1-exp(-theta)-eta+eta*exp(-theta+theta*F(x, ...)))
	return(cdf)

}


qgepg<-function(p, spec, theta=1, eta=0.5, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq((1/theta)*log(((1-exp(-theta)-eta)*p[p>=0&p<=1]+exp(-theta))/((1-eta*p[p>=0&p<=1])*exp(-theta))), ...)
	return(qf)

}


rgepg<-function(n, spec, theta=1, eta=0.5, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qgepg(u, spec, theta=theta, eta=eta, ...)
	return(sf)
}



deepg<-function(x, spec, lambda=1, a=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	pdf<-x
        pdf[log==FALSE]<-a*lambda*f(x, ...)*(F(x, ...))**(a-1)*exp(-lambda*(F(x, ...))**a)/(1-exp(-lambda))
        pdf[log==TRUE]<-log(lambda)+log(a)+fL(x, ...)+(a-1)*FL(x, ...)-lambda*(F(x, ...))**a-log(1-exp(-lambda))
	return(pdf)

}


peepg<-function(x, spec, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-(1-exp(-lambda*(F(x, ...))**a))/(1-exp(-lambda))
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(1-exp(-lambda*(F(x, ...))**a))-log(1-exp(-lambda))
        cdf[log.p==FALSE&lower.tail==FALSE]<-(exp(-lambda*(F(x, ...))**a)-exp(-lambda))/(1-exp(-lambda))
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(exp(-lambda*(F(x, ...))**a)-exp(-lambda))-log(1-exp(-lambda))
	return(cdf)

}


qeepg<-function(p, spec, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq((-(1/lambda)*log(1-p[p>=0&p<=1]*(1-exp(-lambda))))**(1/a), ...)
	return(qf)

}


reepg<-function(n, spec, lambda=1, a=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qeepg(u, spec, lambda=lambda, a=a, ...)
	return(sf)
}



dtessg<-function(x, spec, lambda=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-lambda*f(x, ...)*exp(-lambda*F(x, ...))/(1-exp(-lambda))
        pdf[log==TRUE]<-log(lambda)+fL(x, ...)-lambda*F(x, ...)-log(1-exp(-lambda))
	return(pdf)

}


ptessg<-function(x, spec, lambda=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-(1-exp(-lambda*F(x, ...)))/(1-exp(-lambda))
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(1-exp(-lambda*F(x, ...)))-log(1-exp(-lambda))
        cdf[log.p==FALSE&lower.tail==FALSE]<-(exp(-lambda*F(x, ...))-exp(-lambda))/(1-exp(-lambda))
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(exp(-lambda*F(x, ...))-exp(-lambda))-log(1-exp(-lambda))
	return(cdf)

}


qtessg<-function(p, spec, lambda=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq(-(1/lambda)*log(1-p[p>=0&p<=1]*(1-exp(-lambda))), ...)
	return(qf)

}


rtessg<-function(n, spec, lambda=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qtessg(u, spec, lambda=lambda, ...)
	return(sf)
}



dmog<-function(x, spec, beta=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-beta*f(x, ...)/((beta+(1-beta)*F(x, ...))**2)
        pdf[log==TRUE]<-log(beta)+fL(x, ...)-2*log(beta+(1-beta)*F(x, ...))
	return(pdf)

}

pmog<-function(x, spec, beta=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-F(x, ...)/(beta+(1-beta)*F(x, ...))
        cdf[log.p==TRUE&lower.tail==TRUE]<-FL(x, ...)-log(beta+(1-beta)*F(x, ...))
        cdf[log.p==FALSE&lower.tail==FALSE]<-beta*(1-F(x, ...))/(beta+(1-beta)*F(x, ...))
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(beta)+log(1-F(x, ...))-log(beta+(1-beta)*F(x, ...))
	return(cdf)

}


qmog<-function(p, spec, beta=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq(beta*p[p>=0&p<=1]/(1-(1-beta)*p[p>=0&p<=1]))
	return(qf)

}

rmog<-function(n, spec, beta=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qmog(u, spec, beta=beta, ...)
	return(sf)
}




dexpg<-function(x, spec, a=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	pdf<-x
        pdf[log==FALSE]<-a*f(x, ...)*(F(x, ...))**(a-1)
        pdf[log==TRUE]<-log(a)+fL(x, ...)+(a-1)*FL(x, ...)
	return(pdf)

}


pexpg<-function(x, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-(F(x, ...))**a
        cdf[log.p==TRUE&lower.tail==TRUE]<-a*FL(x, ...)
        cdf[log.p==FALSE&lower.tail==FALSE]<-1-(F(x, ...))**a
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(1-(F(x, ...))**a)
	return(cdf)

}


qexpg<-function(p, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq((p[p>=0&p<=1])**(1/a), ...)
	return(qf)

}


rexpg<-function(n, spec, a=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qexpg(u, spec, a=a, ...)
	return(sf)
}




dexpkumg<-function(x, spec, a=1, b=1, c=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	pdf<-x
        pdf[log==FALSE]<-a*b*c*f(x, ...)*(F(x, ...))**(a-1)*(1-(F(x, ...))**a)**(b-1)*(1-(1-(F(x, ...))**a)**b)**(c-1)
        pdf[log==TRUE]<-log(a)+log(b)+log(c)+fL(x, ...)+(a-1)*FL(x, ...)+(b-1)*log(1-(F(x, ...))**a)+(c-1)*log(1-(1-(F(x, ...))**a)**b)
	return(pdf)

}


pexpkumg<-function(x, spec, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-(1-(1-(F(x, ...))**a)**b)**c
        cdf[log.p==TRUE&lower.tail==TRUE]<-c*log(1-(1-(F(x, ...))**a)**b)
        cdf[log.p==FALSE&lower.tail==FALSE]<-1-(1-(1-(F(x, ...))**a)**b)**c
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(1-(1-(1-(F(x, ...))**a)**b)**c)
	return(cdf)

}


qexpkumg<-function(p, spec, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq((1-(1-(p[p>=0&p<=1])**(1/c))**(1/b))**(1/a), ...)
	return(qf)

}


rexpkumg<-function(n, spec, a=1, b=1, c=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qexpkumg(u, spec, a=a, b=b, c=c, ...)
	return(sf)
}


deg<-function(x, spec, a=1, b=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-a*b*f(x, ...)*(1-F(x, ...))**(a-1)*(1-(1-F(x, ...))**a)**(b-1)
        pdf[log==TRUE]<-log(a)+log(b)+fL(x, ...)+(a-1)*log(1-F(x, ...))+(b-1)*log(1-(1-F(x, ...))**a)
	return(pdf)

}


peg<-function(x, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-(1-(1-F(x, ...))**a)**b
        cdf[log.p==TRUE&lower.tail==TRUE]<-b*log(1-(1-F(x, ...))**a)
        cdf[log.p==FALSE&lower.tail==FALSE]<-1-(1-(1-F(x, ...))**a)**b
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(1-(1-(1-F(x, ...))**a)**b)
	return(cdf)

}


qeg<-function(p, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq(1-(1-(p[p>=0&p<=1])**(1/b))**(1/a), ...)
	return(qf)

}

reg<-function(n, spec, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qeg(u, spec, a=a, b=b, ...)
	return(sf)
}


dkumg<-function(x, spec, a=1, b=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	pdf<-x
        pdf[log==FALSE]<-a*b*f(x, ...)*(F(x, ...))**(a-1)*(1-F(x, ...))**(b-1)
        pdf[log==TRUE]<-log(a)+log(b)+fL(x, ...)+(a-1)*FL(x, ...)+(b-1)*log(1-F(x, ...))
	return(pdf)

}



pkumg<-function(x, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-1-(1-(F(x, ...))**a)**b
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(1-(1-(F(x, ...))**a)**b)
        cdf[log.p==FALSE&lower.tail==FALSE]<-(1-(F(x, ...))**a)**b
        cdf[log.p==TRUE&lower.tail==FALSE]<-b*log(1-(F(x, ...))**a)
	return(cdf)

}


qkumg<-function(p, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq((1-(1-p[p>=0&p<=1])**(1/b))**(1/a), ...)
	return(qf)

}


rkumg<-function(n, spec, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qkumg(u, spec, a=a, b=b, ...)
	return(sf)
}



dbeg<-function(x, spec, alpha=1, a=1, b=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-f(x, ...)*exp(-alpha*F(x, ...))*dbeta(exp(-alpha*F(x, ...)),shape1=a,shape2=b)
        pdf[log==TRUE]<-fL(x, ...)-alpha*F(x, ...)+dbeta(exp(-alpha*F(x, ...)),shape1=a,shape2=b,log=TRUE)
	return(pdf)

}


pbeg<-function(x, spec, alpha=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-pbeta(1-exp(-alpha*F(x, ...)),shape1=a,shape2=b)
        cdf[log.p==TRUE&lower.tail==TRUE]<-pbeta(1-exp(-alpha*F(x, ...)),shape1=a,shape2=b,log.p=TRUE)
        cdf[log.p==FALSE&lower.tail==FALSE]<-pbeta(1-exp(-alpha*F(x, ...)),shape1=a,shape2=b,lower.tail==FALSE)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pbeta(1-exp(-alpha*F(x, ...)),shape1=a,shape2=b,log.p=TRUE,lower.tail==FALSE)
	return(cdf)

}


qbeg<-function(p, spec, alpha=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1-exp(-alpha)]<-Fq(-(1/alpha)*log(1-qbeta(p[p>=0&p<=1-exp(-alpha)],shape1=a,shape2=b)), ...)
	return(qf)

}


rbeg<-function(n, spec, alpha=1, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1-exp(-alpha))
        sf<-qbeg(u, spec, alpha=alpha, a=a, b=b, ...)
	return(sf)
}




dgbg<-function(x, spec, a=1, b=1, c=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	pdf<-x
        pdf[log==FALSE]<-c*f(x, ...)*(F(x, ...))**(c-1)*dbeta((F(x, ...))**c,shape1=a,shape2=b)
        pdf[log==TRUE]<-log(c)+fL(x, ...)+(c-1)*FL(x, ...)+dbeta((F(x, ...))**c,shape1=a,shape2=b,log=TRUE)
	return(pdf)

}


pgbg<-function(x, spec, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-pbeta((F(x, ...))**c,shape1=a,shape2=b)
        cdf[log.p==TRUE&lower.tail==TRUE]<-pbeta((F(x, ...))**c,shape1=a,shape2=b,log.p=TRUE)
        cdf[log.p==FALSE&lower.tail==FALSE]<-pbeta((F(x, ...))**c,shape1=a,shape2=b,lower.tail==FALSE)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pbeta((F(x, ...))**c,shape1=a,shape2=b,log.p=TRUE,lower.tail==FALSE)
	return(cdf)

}


qgbg<-function(p, spec, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq((qbeta(p[p>=0&p<=1],shape1=a,shape2=b))**(1/c), ...)
	return(qf)

}


rgbg<-function(n, spec, a=1, b=1, c=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qgbg(u, spec, a=a, b=b, c=c, ...)
	return(sf)
}

dmbetag<-function(x, spec, beta=1, a=1, b=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-beta**a*f(x, ...)*dbeta(F(x, ...),shape1=a,shape2=b)/(1-(1-beta)*F(x, ...))**(a+b)
        pdf[log==TRUE]<-a*log(beta)+fL(x, ...)+dbeta(F(x, ...),shape1=a,shape2=b,log=TRUE)-(a+b)*log(1-(1-beta)*F(x, ...))
	return(pdf)

}




pmbetag<-function(x, spec, beta=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-pbeta(beta*F(x, ...)/(1-(1-beta)*F(x, ...)),shape1=a,shape2=b)
        cdf[log.p==TRUE&lower.tail==TRUE]<-pbeta(beta*F(x, ...)/(1-(1-beta)*F(x, ...)),shape1=a,shape2=b,log.p=TRUE)
        cdf[log.p==FALSE&lower.tail==FALSE]<-pbeta(beta*F(x, ...)/(1-(1-beta)*F(x, ...)),shape1=a,shape2=b,lower.tail==FALSE)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pbeta(beta*F(x, ...)/(1-(1-beta)*F(x, ...)),shape1=a,shape2=b,log.p=TRUE,lower.tail==FALSE)
	return(cdf)

}



qmbetag<-function(p, spec, beta=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq(qbeta(p[p>=0&p<=1],shape1=a,shape2=b)/(beta-(beta-1)*qbeta(p[p>=0&p<=1],shape1=a,shape2=b)), ...)
	return(qf)

}


rmbetag<-function(n, spec, beta=1, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qmbetag(u, spec, beta=beta, a=a, b=b, ...)
	return(sf)
}




dbetag<-function(x, spec, a=1, b=1, log=FALSE, ...)
{
        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-f(x, ...)*dbeta(F(x, ...),shape1=a,shape2=b)
        pdf[log==TRUE]<-fL(x, ...)+dbeta(F(x, ...),shape1=a,shape2=b,log=TRUE)
	return(pdf)

}


pbetag<-function(x, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-pbeta(F(x, ...),shape1=a,shape2=b)
        cdf[log.p==TRUE&lower.tail==TRUE]<-pbeta(F(x, ...),shape1=a,shape2=b,log.p=TRUE)
        cdf[log.p==FALSE&lower.tail==FALSE]<-pbeta(F(x, ...),shape1=a,shape2=b,lower.tail==FALSE)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pbeta(F(x, ...),shape1=a,shape2=b,log.p=TRUE,lower.tail==FALSE)
	return(cdf)

}

qbetag<-function(p, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq(qbeta(p[p>=0&p<=1],shape1=a,shape2=b), ...)
	return(qf)

}


rbetag<-function(n, spec, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qbetag(u, spec, a=a, b=b, ...)
	return(sf)
}




dloggammag1<-function(x, spec, a=1, b=1, log=FALSE, ...)
{

        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-b**a*f(x, ...)*(-log(1-F(x, ...)))**(a-1)*(1-F(x, ...))**(b-1)/gamma(a)
        pdf[log==TRUE]<-a*log(b)+fL(x, ...)+(a-1)*log(-log(1-F(x, ...)))+(b-1)*log(1-F(x, ...))-lgamma(a)
	return(pdf)

}


ploggammag1<-function(x, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{

        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-1-pgamma(-b*log(1-F(x, ...)),shape=a)
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(1-pgamma(-b*log(1-F(x, ...)),shape=a))
        cdf[log.p==FALSE&lower.tail==FALSE]<-pgamma(-b*log(1-F(x, ...)),shape=a)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pgamma(-b*log(1-F(x, ...)),shape=a,log.p=TRUE)
	return(cdf)

}


qloggammag1<-function(p, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq(1-exp(-(1/b)*qgamma(1-p[p>=0&p<=1],shape=a)), ...)
	return(qf)

}


rloggammag1<-function(n, spec, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qloggammag1(u, spec, a=a, b=b, ...)
	return(sf)
}



dloggammag2<-function(x, spec, a=1, b=1, log=FALSE, ...)
{

        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}
        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	pdf<-x
        pdf[log==FALSE]<-b**a*f(x, ...)*(-FL(x, ...))**(a-1)*(F(x, ...))**(b-1)/gamma(a)
        pdf[log==TRUE]<-a*log(b)+fL(x, ...)+(a-1)*log(-FL(x, ...))+(b-1)*FL(x, ...)-lgamma(a)
	return(pdf)

}


ploggammag2<-function(x, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{

        FL<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, log.p=TRUE, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-pgamma(-b*FL(x, ...),shape=a)
        cdf[log.p==TRUE&lower.tail==TRUE]<-pgamma(-b*FL(x, ...),shape=a,log.p=TRUE)
        cdf[log.p==FALSE&lower.tail==FALSE]<-1-pgamma(-b*FL(x, ...),shape=a)
        cdf[log.p==TRUE&lower.tail==FALSE]<-log(1-pgamma(-b*FL(x, ...),shape=a))
	return(cdf)

}


qloggammag2<-function(p, spec, a=1, b=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq(exp(-(1/b)*qgamma(p[p>=0&p<=1],shape=a)), ...)
	return(qf)

}


rloggammag2<-function(n, spec, a=1, b=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qloggammag2(u, spec, a=a, b=b, ...)
	return(sf)
}




dgammag1<-function(x, spec, a=1, log=FALSE, ...)
{

        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-f(x, ...)*(-log(1-F(x, ...)))**(a-1)/gamma(a)
        pdf[log==TRUE]<-fL(x, ...)+(a-1)*log(-log(1-F(x, ...)))-lgamma(a)
	return(pdf)

}


pgammag1<-function(x, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{

        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-pgamma(-log(1-F(x, ...)),shape=a)
        cdf[log.p==TRUE&lower.tail==TRUE]<-pgamma(-log(1-F(x, ...)),shape=a,log.p=TRUE)
        cdf[log.p==FALSE&lower.tail==FALSE]<-pgamma(-log(1-F(x, ...)),shape=a,lower.tail==FALSE)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pgamma(-log(1-F(x, ...)),shape=a,log.p=TRUE,lower.tail==FALSE)
	return(cdf)

}

qgammag1<-function(p, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq(1-exp(-qgamma(p[p>=0&p<=1],shape=a)), ...)
	return(qf)

}


rgammag1<-function(n, spec, a=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qgammag1(u, spec, a=a, ...)
	return(sf)
}




dgammag2<-function(x, spec, a=1, log=FALSE, ...)
{

        f<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, ...))}
        fL<-function (z, ...) {do.call(paste("d",spec,sep=""),list(z, log=TRUE, ...))}
        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	pdf<-x
        pdf[log==FALSE]<-f(x, ...)*(-log(F(x, ...)))**(a-1)/gamma(a)
        pdf[log==TRUE]<-fL(x, ...)+(a-1)*log(-log(F(x, ...)))-lgamma(a)
	return(pdf)

}

pgammag2<-function(x, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{

        F<-function (z, ...) {do.call(paste("p",spec,sep=""),list(z, ...))}

	cdf<-x
        cdf[log.p==FALSE&lower.tail==TRUE]<-1-pgamma(-log(F(x, ...)),shape=a)
        cdf[log.p==TRUE&lower.tail==TRUE]<-log(1-pgamma(-log(F(x, ...)),shape=a))
        cdf[log.p==FALSE&lower.tail==FALSE]<-pgamma(-log(F(x, ...)),shape=a)
        cdf[log.p==TRUE&lower.tail==FALSE]<-pgamma(-log(F(x, ...)),shape=a,log.p=TRUE)
	return(cdf)

}

qgammag2<-function(p, spec, a=1, log.p=FALSE, lower.tail=TRUE, ...)
{
        if (log.p==TRUE) p<-exp(p)
        if (lower.tail==FALSE) p<-1-p
        Fq<-function (z, ...) {do.call(paste("q",spec,sep=""),list(z, ...))}

	qf<-rep(NaN,length(p))
        qf[p>=0&p<=1]<-Fq(exp(-qgamma(1-p[p>=0&p<=1],shape=a)), ...)
	return(qf)

}

rgammag2<-function(n, spec, a=1, ...)
{	
	u<-runif(n,min=0,max=1)
        sf<-qgammag2(u, spec, a=a, ...)
	return(sf)
}

mbeg<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[4]; dexp(x,rate)}
cum=function(par,x){rate=par[4]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[4]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[4]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[4]; b=par[5]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[4]; rate=par[5]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[4]; rate=par[5]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[4]; dchisq(x,df)}
cum=function(par,x){df=par[4]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[4]; b=par[5]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[4]; b=par[5]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*x)^(-b)}
}


if(g=="log-logistic"){
den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[4]; b=par[5]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[4]; b=par[5]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[4]; b=par[5]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp(x^a)-1))}
}




pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
lambda=par[3]
lambda*f0*beta(a,b)^(-1)*(1-exp(-lambda*c0))^(a-1)*exp(-lambda*b*c0)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
lambda=par[3]
pbeta(1-exp(-lambda*c0),shape1=a,shape2=b)
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}

mbetaexpg<-function(g, data, starts, method="BFGS"){


if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[4]; dexp(x,rate)}
cum=function(par,x){rate=par[4]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[4]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[4]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[4]; b=par[5]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[4]; rate=par[5]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[4]; rate=par[5]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[4]; dchisq(x,df)}
cum=function(par,x){df=par[4]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[4]; b=par[5]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[4]; b=par[5]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[4]; b=par[5]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[4]; b=par[5]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[4]; b=par[5]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp(x^a)-1))}
}


pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
lambda=par[3]
lambda*f0*beta(a,b)^(-1)*(1-c0)^(lambda*b-1)*(1-(1-c0)^lambda)^(a-1)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
lambda=par[3]
1-pbeta((1-c0)**lambda,shape1=lambda*(b-1)+1,shape2=a)
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}


mbetag<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[3]; dexp(x,rate)}
cum=function(par,x){rate=par[3]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[3]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[3]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[3]; b=par[4]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[3]; rate=par[4]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[3]; rate=par[4]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[3]; dchisq(x,df)}
cum=function(par,x){df=par[3]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[3]; b=par[4]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[3]; b=par[4]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*x)^(-b)}
}


if(g=="log-logistic"){
den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[3]; b=par[4]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[3]; b=par[4]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[3]; b=par[4]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp(x^a)-1))}
}

                                           
pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
f0*beta(a,b)^(-1)*c0^(a-1)*(1-c0)^(b-1)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
pbeta(c0,shape1=a,shape2=b)
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}

meepg<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[3]; dexp(x,rate)}
cum=function(par,x){rate=par[3]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[3]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[3]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[3]; b=par[4]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[3]; rate=par[4]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[3]; rate=par[4]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[3]; dchisq(x,df)}
cum=function(par,x){df=par[3]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[3]; b=par[4]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[3]; b=par[4]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[3]; b=par[4]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[3]; b=par[4]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[3]; b=par[4]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp(x^a)-1))}
}


pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
a*b*(1-exp(-b))^(-1)*f0*c0^(a-1)*exp(-b*c0^a)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
(1-exp(-b))^(-1)*(1-exp(-b*c0^a))
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}


meg<-function(g, data, starts, method="BFGS"){


if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[3]; dexp(x,rate)}
cum=function(par,x){rate=par[3]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[3]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[3]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[3]; b=par[4]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[3]; rate=par[4]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[3]; rate=par[4]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[3]; dchisq(x,df)}
cum=function(par,x){df=par[3]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[3]; b=par[4]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[3]; b=par[4]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[3]; b=par[4]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[3]; b=par[4]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[3]; b=par[4]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp(x^a)-1))}
}



pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
r=par[1]
u=par[2]
r*u*f0*((1-c0)^(r-1))*(1-(1-c0)^r)^(u-1)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
r=par[1]
u=par[2]
(1-(1-c0)^r)^u
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else{"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}


mexpg<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[2]; dexp(x,rate)}
cum=function(par,x){rate=par[2]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[2]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[2]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[2]; b=par[3]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[2]; rate=par[3]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[2]; rate=par[3]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[2]; dchisq(x,df)}
cum=function(par,x){df=par[2]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[2]; b=par[3]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[2]; b=par[3]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[2]; b=par[3]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[2]; b=par[3]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[2]; b=par[3]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp(x^a)-1))}
}




pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
a*f0*c0^(a-1)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
c0^a
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}

mexpkumg<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[4]; dexp(x,rate)}
cum=function(par,x){rate=par[4]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[4]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[4]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[4]; b=par[5]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[4]; rate=par[5]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[4]; rate=par[5]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[4]; dchisq(x,df)}
cum=function(par,x){df=par[4]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[4]; b=par[5]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[4]; b=par[5]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[4]; b=par[5]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[4]; b=par[5]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[4]; b=par[5]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp(x^a)-1))}
}


pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
d=par[3]
a*b*d*f0*c0^(a-1)*(1-c0^a)^(b-1)*(1-(1-c0^a)^b)^(d-1)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
d=par[3]
(1-(1-c0^a)^b)^d
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}


mgammag1<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[2]; dexp(x,rate)}
cum=function(par,x){rate=par[2]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[2]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[2]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[2]; b=par[3]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[2]; rate=par[3]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[2]; rate=par[3]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[2]; dchisq(x,df)}
cum=function(par,x){df=par[2]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[2]; b=par[3]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[2]; b=par[3]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[2]; b=par[3]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[2]; b=par[3]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[2]; b=par[3]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp(x^a)-1))}
}

pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
gamma(a)^(-1)*f0*(-log(1-c0))^(a-1)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
pgamma(-log(1-c0),shape=a)
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}


mgammag2<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[2]; dexp(x,rate)}
cum=function(par,x){rate=par[2]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[2]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[2]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[2]; b=par[3]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[2]; rate=par[3]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[2]; rate=par[3]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[2]; dchisq(x,df)}
cum=function(par,x){df=par[2]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[2]; b=par[3]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[2]; b=par[3]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[2]; b=par[3]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[2]; b=par[3]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[2]; b=par[3]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp(x^a)-1))}
}


pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
gamma(a)^(-1)*f0*(-log(c0))^(a-1)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
pgamma(-log(c0),shape=a)
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}


mgammag<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[2]; dexp(x,rate)}
cum=function(par,x){rate=par[2]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[2]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[2]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[2]; b=par[3]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[2]; rate=par[3]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[2]; rate=par[3]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[2]; dchisq(x,df)}
cum=function(par,x){df=par[2]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[2]; b=par[3]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[2]; b=par[3]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[2]; b=par[3]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[2]; b=par[3]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[2]; b=par[3]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp(x^a)-1))}
}



pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
gamma(a)^(-1)*f0*(1-c0)^(-2)*(c0/(1-c0))^(a-1)*exp(-c0/(1-c0))
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
pgamma(c0/(1-c0),shape=a)
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}


mgbg<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[4]; dexp(x,rate)}
cum=function(par,x){rate=par[4]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[4]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[4]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[4]; b=par[5]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[4]; rate=par[5]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[4]; rate=par[5]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[4]; dchisq(x,df)}
cum=function(par,x){df=par[4]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[4]; b=par[5]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[4]; b=par[5]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[4]; b=par[5]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[4]; b=par[5]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[4]; b=par[5]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp(x^a)-1))}
}



pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
d=par[3]
d*beta(a,b)^(-1)*f0*c0^(a*d-1)*(1-c0^d)^(b-1)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
d=par[3]
pbeta(c0^d,shape1=a,shape2=b)
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}

mgepg<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[3]; dexp(x,rate)}
cum=function(par,x){rate=par[3]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[3]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[3]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[3]; b=par[4]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[3]; rate=par[4]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[3]; rate=par[4]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[3]; dchisq(x,df)}
cum=function(par,x){df=par[3]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[3]; b=par[4]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[3]; b=par[4]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[3]; b=par[4]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[3]; b=par[4]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[3]; b=par[4]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp(x^a)-1))}
}




pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
theta=par[1]
eta=par[2]
theta*(1-eta)*(1-exp(-theta))*f0*exp(-theta+theta*c0)/(1-exp(-theta)-eta+eta*exp(-theta+theta*c0))**2
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
theta=par[1]
eta=par[2]
(exp(-theta+theta*c0)-exp(-theta))/(1-exp(-theta)-eta+eta*exp(-theta+theta*c0))
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}


mkumg<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[3]; dexp(x,rate)}
cum=function(par,x){rate=par[3]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[3]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[3]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[3]; b=par[4]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[3]; rate=par[4]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[3]; rate=par[4]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[3]; dchisq(x,df)}
cum=function(par,x){df=par[3]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[3]; b=par[4]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[3]; b=par[4]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[3]; b=par[4]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[3]; b=par[4]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[3]; b=par[4]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp(x^a)-1))}
}





pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
r=par[1]
u=par[2]
r*u*f0*(c0^(r-1))*(1-c0^r)^(u-1)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
r=par[1]
u=par[2]
1-(1-c0^r)^u
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else{"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}


mloggammag1<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[3]; dexp(x,rate)}
cum=function(par,x){rate=par[3]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[3]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[3]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[3]; b=par[4]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[3]; rate=par[4]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[3]; rate=par[4]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[3]; dchisq(x,df)}
cum=function(par,x){df=par[3]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[3]; b=par[4]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[3]; b=par[4]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[3]; b=par[4]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[3]; b=par[4]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[3]; b=par[4]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp(x^a)-1))}
}



pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
b^a*gamma(a)^(-1)*f0*(-log(1-c0))^(a-1)*(1-c0)^(b-1)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
1-pgamma(-b*log(1-c0),shape=a)
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}


mloggammag2<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[3]; dexp(x,rate)}
cum=function(par,x){rate=par[3]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[3]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[3]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[3]; b=par[4]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[3]; rate=par[4]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[3]; rate=par[4]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[3]; dchisq(x,df)}
cum=function(par,x){df=par[3]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[3]; b=par[4]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[3]; b=par[4]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[3]; b=par[4]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[3]; b=par[4]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[3]; b=par[4]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp(x^a)-1))}
}



pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
b^a*gamma(a)^(-1)*f0*(-log(c0))^(a-1)*(c0)^(b-1)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
1-pgamma(-b*log(c0),shape=a)
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}

mmbetag<-function(g, data, starts, method="BFGS"){


if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[4]; dexp(x,rate)}
cum=function(par,x){rate=par[4]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[4]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[4]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[4]; b=par[5]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[4]; scale=par[5]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[4]; scale=par[5]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[4]; rate=par[5]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[4]; rate=par[5]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[4]; sdlog=par[5]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[4]; sdlog=par[5]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[4]; dchisq(x,df)}
cum=function(par,x){df=par[4]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[4]; df2=par[5]; df(x,df1,df2)}
cum=function(par,x){df1=par[4]; df2=par[5]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[4]; d=par[5]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[4]; d=par[5]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[4]; b=par[5]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[4]; b=par[5]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[4]; b=par[5]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[4]; b=par[5]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[4]; b=par[5]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[4]; b=par[5]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[4]; b=par[5]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[4]; b=par[5]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[4]; b=par[5]; 1-exp(-b*(exp(x^a)-1))}
}






pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
d=par[3]
d^a*f0*beta(a,b)^(-1)*c0^(a-1)*(1-c0)^(b-1)*(1-(1-d)*c0)^(-a-b)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
d=par[3]
pbeta(d*c0/(1-(1-d)*c0),shape1=a,shape2=b)
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}


mmog<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[2]; dexp(x,rate)}
cum=function(par,x){rate=par[2]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[2]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[2]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[2]; b=par[3]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[2]; rate=par[3]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[2]; rate=par[3]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[2]; dchisq(x,df)}
cum=function(par,x){df=par[2]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[2]; b=par[3]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[2]; b=par[3]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[2]; b=par[3]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[2]; b=par[3]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[2]; b=par[3]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp(x^a)-1))}
}


pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
r=par[1]
r*f0/((1-(1-r)*(1-c0))^2)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
r=par[1]
1-r*(1-c0)/(1-(1-r)*(1-c0))
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}

mtessg<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[2]; dexp(x,rate)}
cum=function(par,x){rate=par[2]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[2]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[2]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[2]; b=par[3]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[2]; scale=par[3]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[2]; scale=par[3]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[2]; rate=par[3]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[2]; rate=par[3]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[2]; sdlog=par[3]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[2]; sdlog=par[3]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[2]; dchisq(x,df)}
cum=function(par,x){df=par[2]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[2]; df2=par[3]; df(x,df1,df2)}
cum=function(par,x){df1=par[2]; df2=par[3]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[2]; d=par[3]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[2]; d=par[3]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[2]; b=par[3]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[2]; b=par[3]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[2]; b=par[3]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[2]; b=par[3]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[2]; b=par[3]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[2]; b=par[3]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[2]; b=par[3]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[2]; b=par[3]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[2]; b=par[3]; 1-exp(-b*(exp(x^a)-1))}
}


pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
a*(1-exp(-a))^(-1)*f0*exp(-a*c0)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
(1-exp(-a*c0))*(1-exp(-a))^(-1)
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else {"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}



mweibullg<-function(g, data, starts, method="BFGS"){

if(g!="exp" & g!="rayleigh" & g!="weibull" & g!="gompertz" & g!="gamma" 
 & g!="lnorm" & g!="chisq" & g!="f" & g!="burrxii" & g!="frechet" 
 & g!="lomax" & g!="log-logistic" & g!="lfr" & g!="chen")
 { stop ("Baseline distribution not implemented or misspelled. Please check the manual for guidelines.") }

if(g=="exp"){
den=function(par,x){rate=par[3]; dexp(x,rate)}
cum=function(par,x){rate=par[3]; pexp(x,rate)}
}

if(g=="rayleigh"){
den=function(par,x){a=par[3]; 2*x*a*exp(-a*x^2)}
cum=function(par,x){a=par[3]; 1-exp(-a*x^2)}
}

if(g=="gompertz"){
den=function(par,x){a=par[3]; b=par[4]; b*exp(a*x)*exp(-(exp(a*x)-1)*b/a)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-(exp(a*x)-1)*b/a)}
}

if(g=="weibull"){
den=function(par,x){shape=par[3]; scale=par[4]; dweibull(x,shape,scale)}
cum=function(par,x){shape=par[3]; scale=par[4]; pweibull(x,shape,scale)}
}

if(g=="gamma"){
den=function(par,x){shape=par[3]; rate=par[4]; dgamma(x,shape,rate)}
cum=function(par,x){shape=par[3]; rate=par[4]; pgamma(x,shape,rate)}
}

if(g=="lnorm"){
den=function(par,x){meanlog=par[3]; sdlog=par[4]; dlnorm(x,meanlog,sdlog)}
cum=function(par,x){meanlog=par[3]; sdlog=par[4]; plnorm(x,meanlog,sdlog)}
}

if(g=="chisq"){
den=function(par,x){df=par[3]; dchisq(x,df)}
cum=function(par,x){df=par[3]; pchisq(x,df)}
}

if(g=="f"){
den=function(par,x){df1=par[3]; df2=par[4]; df(x,df1,df2)}
cum=function(par,x){df1=par[3]; df2=par[4]; pf(x,df1,df2)}
}

if(g=="burrxii"){
den=function(par,x){a=par[3]; d=par[4]; d*a*(1+(x)^d)^(-a-1)*(x^(d-1))}
cum=function(par,x){a=par[3]; d=par[4]; 1-(1+(x)^d)^(-a)}
}

if(g=="frechet"){
den=function(par,x){a=par[3]; b=par[4]; (a*exp(-(x/b)^(-a))*(x/b)^(-a-1))/(b)}
cum=function(par,x){a=par[3]; b=par[4]; exp(-(x/b)^(-a))}
}

if(g=="lomax"){
den=function(par,x){a=par[3]; b=par[4]; (a*b)/((1+a*x)^(b+1))}
cum=function(par,x){a=par[3]; b=par[4]; 1-(1+a*x)^(-b)}
}

if(g=="log-logistic"){
den=function(par,x){a=par[3]; b=par[4]; (a*b^(-a)*x^(a-1))/(((x/b)^a +1)^2)}
cum=function(par,x){a=par[3]; b=par[4]; 1/((x/b)^(-a)+1)}
}

if(g=="lfr"){
den=function(par,x){a=par[3]; b=par[4]; (a+b*x)*exp(-a*x-(b*x^2)/2)}
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-a*x-(b*x^2)/2)}
}

if(g=="chen"){
den=function(par,x){a=par[3]; b=par[4]; a*b*x^(a-1)*exp(x^a)*exp(-b*(exp(x^a)-1))}     
cum=function(par,x){a=par[3]; b=par[4]; 1-exp(-b*(exp(x^a)-1))}
}



pdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
a*b**(-a)*(f0/(1-c0))*(-log(1-c0))**(a-1)*exp(-b**(-a)*(-log(1-c0))**a)
}

cdf0<-function(par,x){
f0=den(par,x)
c0=cum(par,x)
a=par[1]
b=par[2]
1-exp(-b**(-a)*(-log(1-c0))**a)
}

med=suppressWarnings(goodness.fit(pdf=pdf0, cdf=cdf0, starts=starts, data=data, method=method, mle=NULL))

aux=cbind(med$mle,med$Erro,med$mle+qnorm(0.025)*med$Erro,med$mle+qnorm(0.975)*med$Erro)
colnames(aux)=c("MLE","Std. Dev.","Inf. 95% CI","Sup. 95% CI")

aux1=cbind(med$AIC, med$CAIC, med$BIC, med$HQIC, med$W, med$A, med$Value)
colnames(aux1)=c("AIC","CAIC","BIC","HQIC","W","A", "Min(-log(Likelihood))")
rownames(aux1)=c("")

aux2=cbind(med$KS$statistic,med$KS$p.value)
colnames(aux2)=c("KS Statistic","KS p-value")
rownames(aux2)=c("")

aux3=cbind(if(med$Convergence==0){"Algorithm Converged"} else{"Algorithm Not Converged"})
colnames(aux3)=c("")
rownames(aux3)=c("")

list("Estimates"=aux,"Measures"=aux1,"Kolmogorov-Smirnov Test"=aux2,"Convergence Status"=aux3)
}
