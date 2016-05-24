
#exponential

dexponential=function (x, lambda=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=lambda*exp(-lambda*x)
	pdf[log==TRUE]=log(lambda)-lambda*x
	return(pdf)
}


pexponential=function (x, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=1-exp(-lambda*x)
        cdf[log.p==FALSE&lower.tail==FALSE]=exp(-lambda*x)
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-exp(-lambda*x))
        cdf[log.p==TRUE&lower.tail==FALSE]=-lambda*x
	return(cdf)
}


varexponential=function (p, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(-1/lambda)*log(1-p)
	return(var)
}


esexponential=function (p, lambda=1)
{
	f=function (x) {varexponential(x,lambda=lambda)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}


#Kumaraswamy exponential

dkumexp=function (x, lambda=1, a=1, b=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=a*b*lambda*exp(-lambda*x)*(1-exp(-lambda*x))**(a-1)*(1-(1-exp(-lambda*x))**a)**(b-1)
	pdf[log==TRUE]=log(a*b*lambda)-lambda*x+(a-1)*log(1-exp(-lambda*x))+(b-1)*log(1-(1-exp(-lambda*x))**a)
	return(pdf)
}

pkumexp=function (x, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=1-(1-(1-exp(-lambda*x))**a)**b
        cdf[log.p==FALSE&lower.tail==FALSE]=(1-(1-exp(-lambda*x))**a)**b
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-(1-(1-exp(-lambda*x))**a)**b)
        cdf[log.p==TRUE&lower.tail==FALSE]=b*log(1-(1-exp(-lambda*x))**a)
	return(cdf)
}

varkumexp=function (p, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(-1/lambda)*log(1-(1-(1-p)**(1/b))**(1/a))
	return(var)
}

eskumexp=function (p, lambda=1, a=1, b=1)
{
	f=function (x) {varkumexp(x,lambda=lambda,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}


#Exponentiated exponential

dexpexp=function (x, lambda=1, a=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=a*lambda*exp(-lambda*x)*(1-exp(-lambda*x))**(a-1)
	pdf[log==TRUE]=log(a*lambda)-lambda*x+(a-1)*log(1-exp(-lambda*x))
	return(pdf)
}

pexpexp=function (x, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=(1-exp(-lambda*x))**a
        cdf[log.p==FALSE&lower.tail==FALSE]=1-(1-exp(-lambda*x))**a
        cdf[log.p==TRUE&lower.tail==TRUE]=a*log(1-exp(-lambda*x))
        cdf[log.p==TRUE&lower.tail==FALSE]=log(1-(1-exp(-lambda*x))**a)
	return(cdf)
}

varexpexp=function (p, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(-1/lambda)*log(1-p**(1/a))
	return(var)
}

esexpexp=function (p, lambda=1, a=1)
{
	f=function (x) {varexpexp(x,lambda=lambda,a=a)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Inverse exponentiated exponential

dinvexpexp=function (x, lambda=1, a=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=a*lambda*(1/(x*x))*exp(-lambda/x)*(1-exp(-lambda/x))**(a-1)
	pdf[log==TRUE]=log(a*lambda)-2*log(x)-lambda/x+(a-1)*log(1-exp(-lambda/x))
	return(pdf)
}

pinvexpexp=function (x, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=1-(1-exp(-lambda/x))**a
        cdf[log.p==FALSE&lower.tail==FALSE]=(1-exp(-lambda/x))**a
        cdf[log.p==TRUE&lower.tail==TRUE]=a*log(1-(1-exp(-lambda/x))**a)
        cdf[log.p==TRUE&lower.tail==FALSE]=a*log(1-exp(-lambda/x))
	return(cdf)
}

varinvexpexp=function (p, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=lambda/(-log(1-(1-p)**(1/a)))
	return(var)
}

esinvexpexp=function (p, lambda=1, a=1)
{
	f=function (x) {varinvexpexp(x,lambda=lambda,a=a)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Beta exponential

dbetaexp=function (x, lambda=1, a=1, b=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=lambda*exp(-b*lambda*x)*(1-exp(-lambda*x))**(a-1)/beta(a,b)
	pdf[log==TRUE]=log(lambda)-b*lambda*x+(a-1)*log(1-exp(-lambda*x))-lbeta(a,b)
	return(pdf)
}


pbetaexp=function (x, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(1-exp(-lambda*x),shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}


varbetaexp=function (p, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(-1/lambda)*log(1-qbeta(p,shape1=a,shape2=b))
	return(var)
}

esbetaexp=function (p, lambda=1, a=1, b=1)
{
	f=function (x) {varbetaexp(x,lambda=lambda,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Logistic exponential

dlogisexp=function (x, lambda=1, a=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=a*lambda*exp(lambda*x)*(exp(lambda*x)-1)**(a-1)/(1+(exp(lambda*x)-1)**a)**2
        pdf[log==TRUE]=log(a*lambda)+lambda*x+(a-1)*log(exp(lambda*x)-1)-2*log(1+(exp(lambda*x)-1)**a)
	return(pdf)
}

plogisexp=function (x, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=(exp(lambda*x)-1)**a/(1+(exp(lambda*x)-1)**a)
        cdf[log.p==FALSE&lower.tail==FALSE]=1/(1+(exp(lambda*x)-1)**a)
        cdf[log.p==TRUE&lower.tail==TRUE]=a*log(exp(lambda*x)-1)-log(1+(exp(lambda*x)-1)**a)
        cdf[log.p==TRUE&lower.tail==FALSE]=-log(1+(exp(lambda*x)-1)**a)
	return(cdf)
}

varlogisexp=function (p, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(1/lambda)*log(1+(p/(1-p))**(1/a))
	return(var)
}

eslogisexp=function (p, lambda=1, a=1)
{
	f=function (x) {varlogisexp(x,lambda=lambda,a=a)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Exponential extension

dexpext=function (x, lambda=1, a=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=a*lambda*(1+lambda*x)**(a-1)*exp(1-(1+lambda*x)**a)
        pdf[log==TRUE]=log(a*lambda)+(a-1)*log(1+lambda*x)+1-(1+lambda*x)**a
	return(pdf)
}

pexpext=function (x, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=1-exp(1-(1+lambda*x)**a)
        cdf[log.p==FALSE&lower.tail==FALSE]=exp(1-(1+lambda*x)**a)
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-exp(1-(1+lambda*x)**a))
        cdf[log.p==TRUE&lower.tail==FALSE]=1-(1+lambda*x)**a
	return(cdf)
}

varexpext=function (p, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(1/lambda)*((1-log(1-p))**(1/a)-1)
	return(var)
}

esexpext=function (p, lambda=1, a=1)
{
	f=function (x) {varexpext(x,lambda=lambda,a=a)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Marshall Olkin exponential

dmoexp=function (x, lambda=1, a=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=lambda*exp(lambda*x)/(exp(lambda*x)-1+a)**2
        pdf[log==TRUE]=log(lambda)+lambda*x-2*log(exp(lambda*x)-1+a)
	return(pdf)
}

pmoexp=function (x, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=(exp(lambda*x)-2+a)/(exp(lambda*x)-1+a)
        cdf[log.p==FALSE&lower.tail==FALSE]=1/(exp(lambda*x)-1+a)
        cdf[log.p==TRUE&lower.tail==TRUE]=log(exp(lambda*x)-2+a)-log(exp(lambda*x)-1+a)
        cdf[log.p==TRUE&lower.tail==FALSE]=-log(exp(lambda*x)-1+a)
	return(cdf)
}

varmoexp=function (p, lambda=1, a=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(1/lambda)*(log(2-a-(1-a)*p)-log(1-p))
	return(var)
}

esmoexp=function (p, lambda=1, a=1)
{
	f=function (x) {varmoexp(x,lambda=lambda,a=a)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Perks

dperks=function (x, a=1, b=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=a*(a+1)*exp(b*x)/(1+a*exp(b*x))**2
        pdf[log==TRUE]=log(a*(a+1))+b*x-2*log(1+a*exp(b*x))
	return(pdf)
}

pperks=function (x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=a*(exp(b*x)-1)/(1+a*exp(b*x))
        cdf[log.p==FALSE&lower.tail==FALSE]=(1+a)/(1+a*exp(b*x))
        cdf[log.p==TRUE&lower.tail==TRUE]=log(a)+log(exp(b*x)-1)-log(1+a*exp(b*x))
        cdf[log.p==TRUE&lower.tail==FALSE]=log(1+a)-log(1+a*exp(b*x))
	return(cdf)
}


varperks=function (p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(1/b)*(log(a+p)-log(a)-log(1-p))
	return(var)
}


esperks=function (p, a=1, b=1)
{
	f=function (x) {varperks(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}


#Beard

dbeard=function (x, a=1, b=1, rho=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=a*exp(b*x)*(1+a*rho)**(rho**(-1/b))/(1+a*rho*exp(b*x))**(1+rho**(-1/b))
        pdf[log==TRUE]=log(a)+b*x+(rho**(-1/b))*log(1+a*rho)-(1+rho**(-1/b))*log(1+a*rho*exp(b*x))
	return(pdf)
}

pbeard=function (x, a=1, b=1, rho=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=1-((1+a*rho)/(1+a*rho*exp(b*x)))**(rho**(-1/b))
        cdf[log.p==FALSE&lower.tail==FALSE]=((1+a*rho)/(1+a*rho*exp(b*x)))**(rho**(-1/b))
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-((1+a*rho)/(1+a*rho*exp(b*x)))**(rho**(-1/b)))
        cdf[log.p==TRUE&lower.tail==FALSE]=(rho**(-1/b))*(log(1+a*rho)-log(1+a*rho*exp(b*x)))
	return(cdf)
}

varbeard=function (p, a=1, b=1, rho=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(1/b)*(log((1+a*rho)/(1-p)**(rho**(1/b))-1)-log(a*rho))
	return(var)
}


esbeard=function (p, a=1, b=1, rho=1)
{
	f=function (x) {varbeard(x,a=a,b=b,rho=rho)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}


#Gompertz

dgompertz=function (x, b=1, eta=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=b*eta*exp(b*x)*exp(eta-eta*exp(b*x))
        pdf[log==TRUE]=log(b*eta)+b*x+eta-eta*exp(b*x)
	return(pdf)
}


pgompertz=function (x, b=1, eta=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=1-exp(eta-eta*exp(b*x))
        cdf[log.p==FALSE&lower.tail==FALSE]=exp(eta-eta*exp(b*x))
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-exp(eta-eta*exp(b*x)))
        cdf[log.p==TRUE&lower.tail==FALSE]=eta-eta*exp(b*x)
	return(cdf)
}

vargompertz=function (p, b=1, eta=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(1/b)*log(1-(1/eta)*log(1-p))
	return(var)
}

esgompertz=function (p, b=1, eta=1)
{
	f=function (x) {vargompertz(x,b=b,eta=eta)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Beta Gompertz

dbetagompertz=function (x, b=1, c=1, d=1, eta=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=b*eta*exp(d*eta-d*eta*exp(b*x)+b*x)*(1-exp(eta-eta*exp(b*x)))**(c-1)
        pdf[log==TRUE]=log(b*eta)+d*eta-d*eta*exp(b*x)+b*x+(c-1)*log(1-exp(eta-eta*exp(b*x)))
	return(pdf)
}

pbetagompertz=function (x, b=1, c=1, d=1, eta=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(1-exp(eta-eta*exp(b*x)),shape1=c,shape2=d,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varbetagompertz=function (p, b=1, c=1, d=1, eta=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(1/b)*log(1-(1/eta)*log(1-qbeta(p,shape1=c,shape2=d)))
	return(var)
}

esbetagompertz=function (p, b=1, c=1, d=1, eta=1)
{
	f=function (x) {varbetagompertz(x,b=b,c=c,d=d,eta=eta)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Linear failure rate


dlfr=function (x, a=1, b=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=(a+b*x)*exp(-a*x-b*x*x/2)
        pdf[log==TRUE]=log(a+b*x)-a*x-b*x*x/2
	return(pdf)
}

plfr=function (x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=1-exp(-a*x-b*x*x/2)
        cdf[log.p==FALSE&lower.tail==FALSE]=exp(-a*x-b*x*x/2)
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-exp(-a*x-b*x*x/2))
        cdf[log.p==TRUE&lower.tail==FALSE]=-a*x-b*x*x/2
	return(cdf)
}

varlfr=function (p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(1/b)*(sqrt(a*a-2*b*log(1-p))-a)
	return(var)
}

eslfr=function (p, a=1, b=1)
{
	f=function (x) {varlfr(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Pareto

dpareto=function (x, K=1, c=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=c*K**c*x**(-c-1)
        pdf[log==TRUE]=log(c)+c*log(K)-(c+1)*log(x)
	return(pdf)
}

ppareto=function (x, K=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=1-K**c*x**(-c)
        cdf[log.p==FALSE&lower.tail==FALSE]=K**c*x**(-c)
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-K**c*x**(-c))
        cdf[log.p==TRUE&lower.tail==FALSE]=c*(log(K)-log(x))
	return(cdf)
}

varpareto=function (p, K=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=K*(1-p)**(-1/c)
	return(var)
}

espareto=function (p, K=1, c=1)
{
	f=function (x) {varpareto(x,K=K,c=c)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Kumaraswamy Pareto

dkumpareto=function (x, K=1, a=1, b=1, c=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=a*b*c*K**c*x**(-c-1)*(1-K**c*x**(-c))**(a-1)*(1-(1-K**c*x**(-c))**a)**(b-1)
        pdf[log==TRUE]=log(a*b*c)+c*log(K)-(c+1)*log(x)+(a-1)*log(1-K**c*x**(-c))+(b-1)*log(1-(1-K**c*x**(-c))**a)
	return(pdf)
}

pkumpareto=function (x, K=1, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=1-(1-(1-K**c*x**(-c))**a)**b
        cdf[log.p==FALSE&lower.tail==FALSE]=(1-(1-K**c*x**(-c))**a)**b
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-(1-(1-K**c*x**(-c))**a)**b)
        cdf[log.p==TRUE&lower.tail==FALSE]=b*log(1-(1-K**c*x**(-c))**a)
	return(cdf)
}

varkumpareto=function (p, K=1, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=K*(1-(1-(1-p)**(1/b))**(1/a))**(-1/c)
	return(var)
}

eskumpareto=function (p, K=1, a=1, b=1, c=1)
{
	f=function (x) {varkumpareto(x,K=K,a=a,b=b,c=c)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#F

dF=function (x, d1=1, d2=1, log=FALSE)
{
	pdf=df(x,df1=d1,df2=d2,log=log)
	return(pdf)
}

pF=function (x, d1=1, d2=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pf(x,df1=d1,df2=d2,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varF=function (p, d1=1, d2=1, log.p=FALSE, lower.tail=TRUE)
{
	var=qf(p,df1=d1,df2=d2,log.p=log.p,lower.tail=lower.tail)
	return(var)
}

esF=function (p, d1=1, d2=1)
{
	f=function (x) {varF(x,d1=d1,d2=d2)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Generalized Pareto


dgenpareto=function (x, k=1, c=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=(1/k)*(1-c*x/k)**(1/c-1)
        pdf[log==TRUE]=-log(k)+(1/c-1)*log(1-c*x/k)
	return(pdf)
}

pgenpareto=function (x, k=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=1-(1-c*x/k)**(1/c)
        cdf[log.p==FALSE&lower.tail==FALSE]=(1-c*x/k)**(1/c)
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-(1-c*x/k)**(1/c))
        cdf[log.p==TRUE&lower.tail==FALSE]=(1/c)*log(1-c*x/k)
	return(cdf)
}

vargenpareto=function (p, k=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(k/c)*(1-(1-p)**c)
	return(var)
}

esgenpareto=function (p, k=1, c=1)
{
	f=function (x) {vargenpareto(x,k=k,c=c)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Beta Pareto


dbetapareto=function (x, K=1, a=1, c=1, d=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=a*K**(a*d)*x**(-a*d-1)*(1-(K/x)**a)**(c-1)/beta(c,d)
        pdf[log==TRUE]=log(a)+(a*d)*log(K)-(a*d+1)*log(x)+(c-1)*log(1-(K/x)**a)-lbeta(c,d)
	return(pdf)
}

pbetapareto=function (x, K=1, a=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(1-(K/x)**a,shape1=c,shape2=d,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varbetapareto=function (p, K=1, a=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=K*(1-qbeta(p,shape1=c,shape2=d))**(-1/a)
	return(var)
}

esbetapareto=function (p, K=1, a=1, c=1, d=1)
{
	f=function (x) {varbetapareto(x,K=K,a=a,c=c,d=d)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Pareto positive stable


dparetostable=function (x, lambda=1, nu=1, sigma=1, log=FALSE)
{
	pdf=x
	pdf[log==FALSE]=nu*lambda*(1/x)*(log(x/sigma))**(nu-1)*exp(-lambda*(log(x/sigma))**nu)
        pdf[log==TRUE]=log(nu*lambda)-log(x)+(nu-1)*log(log(x/sigma))-lambda*(log(x/sigma))**nu
	return(pdf)
}

pparetostable=function (x, lambda=1, nu=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=1-exp(-lambda*(log(x/sigma))**nu)
        cdf[log.p==FALSE&lower.tail==FALSE]=exp(-lambda*(log(x/sigma))**nu)
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-exp(-lambda*(log(x/sigma))**nu))
        cdf[log.p==TRUE&lower.tail==FALSE]=-lambda*(log(x/sigma))**nu
	return(cdf)
}

varparetostable=function (p, lambda=1, nu=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=sigma*exp(((-1/lambda)*log(1-p))**(1/nu))
	return(var)
}

esparetostable=function (p, lambda=1, nu=1, sigma=1)
{
	f=function (x) {varparetostable(x,lambda=lambda,nu=nu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Gamma


dGamma=function (x, a=1, b=1, log=FALSE)
{
	pdf=dgamma(x,shape=a,rate=b,log=log)
	return(pdf)
}

pGamma=function (x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pgamma(x,shape=a,rate=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varGamma=function (p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	var=qgamma(p,shape=a,rate=b,log.p=log.p,lower.tail=lower.tail)
	return(var)
}

esGamma=function (p, a=1, b=1)
{
	f=function (x) {varGamma(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Kumaraswamy gamma

dkumgamma=function (x, a=1, b=1, c=1, d=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=c*d*dgamma(x,shape=a,rate=b)*(pgamma(x,shape=a,rate=b))**(c-1)*(1-(pgamma(x,shape=a,rate=b))**c)**(d-1)
        pdf[log==TRUE]=log(c*d)+dgamma(x,shape=a,rate=b,log=TRUE)+(c-1)*pgamma(x,shape=a,rate=b,log.p=TRUE)+(d-1)*log(1-(pgamma(x,shape=a,rate=b))**c)
	return(pdf)
}

pkumgamma=function (x, a=1, b=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=1-(1-(pgamma(x,shape=a,rate=b))**c)**d
        cdf[log.p==FALSE&lower.tail==FALSE]=(1-(pgamma(x,shape=a,rate=b))**c)**d
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-(1-(pgamma(x,shape=a,rate=b))**c)**d)
        cdf[log.p==TRUE&lower.tail==FALSE]=d*log(1-(pgamma(x,shape=a,rate=b))**c)
	return(cdf)
}

varkumgamma=function (p, a=1, b=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(1/b)*qgamma((1-(1-p)**(1/d))**(1/c),shape=a)
	return(var)
}

eskumgamma=function (p, a=1, b=1, c=1, d=1)
{
	f=function (x) {varkumgamma(x,a=a,b=b,c=c,d=d)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Nakagami

dnakagami=function (x, m=1, a=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=2*(m/a)**m*x**(2*m-1)*exp(-m*x*x/a)/gamma(m)
        pdf[log==TRUE]=log(2)+m*log(m/a)+(2*m-1)*log(x)-m*x*x/a-lgamma(m)
	return(pdf)
}

pnakagami=function (x, m=1, a=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=pgamma(m*x*x/a,shape=m,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varnakagami=function (p, m=1, a=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=sqrt(a/m)*sqrt(qgamma(p,shape=m))
	return(var)
}

esnakagami=function (p, m=1, a=1)
{
	f=function (x) {varnakagami(x,m=m,a=a)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Reflected gamma

drgamma=function (x, a=1, theta=0, phi=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(1/(2*phi))*dgamma(abs(x-theta)/phi,shape=a)
        pdf[log==TRUE]=dgamma(abs(x-theta)/phi,shape=a,log=TRUE)-log(2*phi)
	return(pdf)
}

prgamma=function (x, a=1, theta=0, phi=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[x<=theta&log.p==FALSE&lower.tail==TRUE]=0.5-0.5*pgamma((theta-x)/phi,shape=a)
        cdf[x>theta&log.p==FALSE&lower.tail==TRUE]=0.5+0.5*pgamma((x-theta)/phi,shape=a)
	cdf[x<=theta&log.p==FALSE&lower.tail==FALSE]=0.5+0.5*pgamma((theta-x)/phi,shape=a)
        cdf[x>theta&log.p==FALSE&lower.tail==FALSE]=0.5-0.5*pgamma((x-theta)/phi,shape=a)
	cdf[x<=theta&log.p==TRUE&lower.tail==TRUE]=log(0.5-0.5*pgamma((theta-x)/phi,shape=a))
        cdf[x>theta&log.p==TRUE&lower.tail==TRUE]=log(0.5+0.5*pgamma((x-theta)/phi,shape=a))
	cdf[x<=theta&log.p==TRUE&lower.tail==FALSE]=log(0.5+0.5*pgamma((theta-x)/phi,shape=a))
        cdf[x>theta&log.p==TRUE&lower.tail==FALSE]=log(0.5-0.5*pgamma((x-theta)/phi,shape=a))
	return(cdf)
}

varrgamma=function (p, a=1, theta=0, phi=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=p
        var[p<=0.5]=theta-phi*qgamma(1-2*p[p<=0.5],shape=a)
        var[p>0.5]=theta+phi*qgamma(2*p[p>0.5]-1,shape=a)
	return(var)
}

esrgamma=function (p, a=1, theta=0, phi=1)
{
	f=function (x) {varrgamma(x,a=a,theta=theta,phi=phi)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}


#Compound Laplace gamma

dclg=function (x, a=1, b=1, theta=0, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=0.5*a*b*(1+b*abs(x-theta))**(-a-1)
        pdf[log==TRUE]=log(a*b)-log(2)-(a+1)*log(1+b*abs(x-theta))
	return(pdf)
}

pclg=function (x, a=1, b=1, theta=0, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[x<=theta&log.p==FALSE&lower.tail==TRUE]=0.5*(1+b*abs(x-theta))**(-a)
        cdf[x>theta&log.p==FALSE&lower.tail==TRUE]=1-0.5*(1+b*abs(x-theta))**(-a)
	cdf[x<=theta&log.p==FALSE&lower.tail==FALSE]=1-0.5*(1+b*abs(x-theta))**(-a)
        cdf[x>theta&log.p==FALSE&lower.tail==FALSE]=0.5*(1+b*abs(x-theta))**(-a)
	cdf[x<=theta&log.p==TRUE&lower.tail==TRUE]=-log(2)-a*log(1+b*abs(x-theta))
        cdf[x>theta&log.p==TRUE&lower.tail==TRUE]=log(1-0.5*(1+b*abs(x-theta))**(-a))
	cdf[x<=theta&log.p==TRUE&lower.tail==FALSE]=log(1-0.5*(1+b*abs(x-theta))**(-a))
        cdf[x>theta&log.p==TRUE&lower.tail==FALSE]=-log(2)-a*log(1+b*abs(x-theta))
	return(cdf)
}

varclg=function (p, a=1, b=1, theta=0, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=p
        var[p<=0.5]=theta-1/b-(1/b)*(2*p[p<=0.5])**(-1/a)
        var[p>0.5]=theta-1/b+(1/b)*(2*(1-p[p>0.5]))**(-1/a)
	return(var)
}

esclg=function (p, a=1, b=1, theta=0)
{
	f=function (x) {varclg(x,a=a,b=b,theta=theta)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Log gamma


dloggamma=function (x, a=1, r=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(1/x)*dgamma(-log(x),shape=r,rate=a)
        pdf[log==TRUE]=dgamma(-log(x),shape=r,rate=a,log=TRUE)-log(x)
	return(pdf)
}

ploggamma=function (x, a=1, r=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
	cdf[log.p==FALSE&lower.tail==TRUE]=1-pgamma(-a*log(x),shape=r)
        cdf[log.p==FALSE&lower.tail==FALSE]=pgamma(-a*log(x),shape=r)
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-pgamma(-a*log(x),shape=r))
        cdf[log.p==TRUE&lower.tail==FALSE]=pgamma(-a*log(x),shape=r,log.p=TRUE)
	return(cdf)
}

varloggamma=function (p, a=1, r=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=exp(-(1/a)*qgamma(1-p,shape=r))
	return(var)
}

esloggamma=function (p, a=1, r=1)
{
	f=function (x) {varloggamma(x,a=a,r=r)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}


#Inverse gamma


dinvgamma=function (x, a=1, b=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(1/x)*dgamma(1/x,shape=a,rate=b)
        pdf[log==TRUE]=dgamma(1/x,shape=a,rate=b,log=TRUE)-log(x)
	return(pdf)
}

pinvgamma=function (x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-pgamma(b/x,shape=a)
        cdf[log.p==FALSE&lower.tail==FALSE]=pgamma(b/x,shape=a)
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-pgamma(b/x,shape=a))
        cdf[log.p==TRUE&lower.tail==FALSE]=pgamma(b/x,shape=a,log.p=TRUE)
	return(cdf)
}

varinvgamma=function (p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=b/qgamma(1-p,shape=a)
	return(var)
}

esinvgamma=function (p, a=1, b=1)
{
	f=function (x) {varinvgamma(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Stacy gamma


dstacygamma=function (x, gamma=1, c=1, theta=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(c*x**(c-1)/theta**c)*dgamma((x/theta)**c,shape=gamma)
        pdf[log==TRUE]=log(c)+(c-1)*log(x)-c*log(theta)+dgamma((x/theta)**c,shape=gamma,log=TRUE)
	return(pdf)
}

pstacygamma=function (x, gamma=1, c=1, theta=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pgamma((x/theta)**c,shape=gamma,log.p=log.p,lower.tail=TRUE)
	return(cdf)
}

varstacygamma=function (p, gamma=1, c=1, theta=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=theta*(qgamma(p,shape=gamma))**(1/c)
	return(var)
}

esstacygamma=function (p, gamma=1, c=1, theta=1)
{
	f=function (x) {varstacygamma(x,gamma=gamma,c=c,theta=theta)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}


#beta

dbetadist=function (x, a=1, b=1, log=FALSE)
{
	pdf=dbeta(x,shape1=a,shape2=b,log=log)
	return(pdf)
}

pbetadist=function (x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(x,shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varbetadist=function (p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	var=qbeta(p,shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(var)
}

esbetadist=function (p, a=1, b=1)
{
	f=function (x) {varbetadist(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}


#Uniform

duniform=function (x, a=0, b=1, log=FALSE)
{
	pdf=dunif(x,min=a,max=b,log=log)
	return(pdf)
}

puniform=function (x, a=0, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=punif(x,min=a,max=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varuniform=function (p, a=0, b=1, log.p=FALSE, lower.tail=TRUE)
{
	var=qunif(p,min=a,max=b,log.p=log.p,lower.tail=lower.tail)
	return(var)
}

esuniform=function (p, a=0, b=1)
{
	f=function (x) {varuniform(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Generalized uniform

dgenunif=function (x, a=0, c=1, h=1, k=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=h*k*c*(x-a)**(c-1)*(1-k*(x-a)**c)**(h-1)
        pdf[log==TRUE]=log(h*k*c)+(c-1)*log(x-a)+(h-1)*log(1-k*(x-a)**c)
	return(pdf)
}

pgenunif=function (x, a=0, c=1, h=1, k=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-(1-k*(x-a)**c)**h
        cdf[log.p==FALSE&lower.tail==FALSE]=(1-k*(x-a)**c)**h
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-(1-k*(x-a)**c)**h)
        cdf[log.p==TRUE&lower.tail==FALSE]=h*log(1-k*(x-a)**c)
	return(cdf)
}

vargenunif=function (p, a=0, c=1, h=1, k=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=a+k**(-1/c)*(1-(1-p)**(1/h))**(1/c)
	return(var)
}

esgenunif=function (p, a=0, c=1, h=1, k=1)
{
	f=function (x) {vargenunif(x,a=a,c=c,h=h,k=k)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#power function I

dpower1=function (x, a=1, log=FALSE)
{
	pdf=dbeta(x,shape1=a,shape2=1,log=log)
	return(pdf)
}

ppower1=function (x, a=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(x,shape1=a,shape2=1,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varpower1=function (p, a=1, log.p=FALSE, lower.tail=TRUE)
{
	var=qbeta(p,shape1=a,shape2=1,log.p=log.p,lower.tail=lower.tail)
	return(var)
}

espower1=function (p, a=1)
{
	f=function (x) {varpower1(x,a=a)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#power function II

dpower2=function (x, b=1, log=FALSE)
{
	pdf=dbeta(x,shape1=1,shape2=b,log=log)
	return(pdf)
}

ppower2=function (x, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(x,shape1=1,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varpower2=function (p, b=1, log.p=FALSE, lower.tail=TRUE)
{
	var=qbeta(p,shape1=1,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(var)
}

espower2=function (p, b=1)
{
	f=function (x) {varpower2(x,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Log beta

dlogbeta=function (x, a=1, b=1, c=1, d=2, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(dbeta(log(x/c)/log(d/c),shape1=a,shape2=b))/(x*log(d/c))
        pdf[log==TRUE]=dbeta(log(x/c)/log(d/c),shape1=a,shape2=b,log=TRUE)-log(x)-log(log(d/c))
	return(pdf)
}

plogbeta=function (x, a=1, b=1, c=1, d=2, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(log(x/c)/log(d/c),shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varlogbeta=function (p, a=1, b=1, c=1, d=2, log.p=FALSE, lower.tail=TRUE)
{
	var=c*(d/c)**(qbeta(p,shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail))
	return(var)
}

eslogbeta=function (p, a=1, b=1, c=1, d=2)
{
	f=function (x) {varlogbeta(x,a=a,b=b,c=c,d=d)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Complementary beta

dcompbeta=function (x, a=1, b=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=beta(a,b)*(qbeta(x,shape1=a,shape2=b))**(1-a)*(1-qbeta(x,shape1=a,shape2=b))**(1-b)
        pdf[log==TRUE]=lbeta(a,b)+(1-a)*log(qbeta(x,shape1=a,shape2=b))+(1-b)*log(1-qbeta(x,shape1=a,shape2=b))
	return(pdf)
}

pcompbeta=function (x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=qbeta(x,shape1=a,shape2=b)
        cdf[log.p==FALSE&lower.tail==FALSE]=1-qbeta(x,shape1=a,shape2=b)
        cdf[log.p==TRUE&lower.tail==TRUE]=log(qbeta(x,shape1=a,shape2=b))
        cdf[log.p==TRUE&lower.tail==FALSE]=log(1-qbeta(x,shape1=a,shape2=b))
	return(cdf)
}

varcompbeta=function (p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=pbeta(p,shape1=a,shape2=b)
	return(var)
}

escompbeta=function (p, a=1, b=1)
{
	f=function (x) {varcompbeta(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Libby Novick beta

dLNbeta=function (x, lambda=1, a=1, b=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=lambda*(1+(lambda-1)*x)**(-2)*dbeta(lambda*x/(1+(lambda-1)*x),shape1=a,shape2=b)
        pdf[log==TRUE]=log(lambda)-2*log(1+(lambda-1)*x)+dbeta(lambda*x/(1+(lambda-1)*x),shape1=a,shape2=b,log=TRUE)
	return(pdf)
}

pLNbeta=function (x, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(lambda*x/(1+(lambda-1)*x),shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}


varLNbeta=function (p, lambda=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=qbeta(p,shape1=a,shape2=b)/(lambda-(lambda-1)*qbeta(p,shape1=a,shape2=b))
	return(var)
}

esLNbeta=function (p, lambda=1, a=1, b=1)
{
	f=function (x) {varLNbeta(x,lambda=lambda,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}


#McDonald Richards beta

dMRbeta=function (x, a=1, b=1, r=1, q=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(r/b)*(x**(r-1)/q**r)*dbeta((1/b)*(x/q)**r,shape1=a,shape2=b)
        pdf[log==TRUE]=log(r)-log(b)-r*log(q)+(r-1)*log(x)+dbeta((1/b)*(x/q)**r,shape1=a,shape2=b,log=TRUE)
	return(pdf)
}

pMRbeta=function (x, a=1, b=1, r=1, q=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta((1/b)*(x/q)**r,shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}


varMRbeta=function (p, a=1, b=1, r=1, q=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=b**(1/r)*q*(qbeta(p,shape1=a,shape2=b))**(1/r)
	return(var)
}

esMRbeta=function (p, a=1, b=1, r=1, q=1)
{
	f=function (x) {varMRbeta(x,a=a,b=b,r=r,q=q)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}


#Generalized beta

dgenbeta=function (x, a=1, b=1, c=0, d=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(1/(d-c))*dbeta((x-c)/(d-c),shape1=a,shape2=b)
        pdf[log==TRUE]=dbeta((x-c)/(d-c),shape1=a,shape2=b,log=TRUE)-log(d-c)
	return(pdf)
}

pgenbeta=function (x, a=1, b=1, c=0, d=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta((x-c)/(d-c),shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}


vargenbeta=function (p, a=1, b=1, c=0, d=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=c+(d-c)*qbeta(p,shape1=a,shape2=b)
	return(var)
}

esgenbeta=function (p, a=1, b=1, c=0, d=1)
{
	f=function (x) {vargenbeta(x,a=a,b=b,c=c,d=d)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}


#Arcsine

darcsine=function (x, a=0, b=1, log=FALSE)
{
	pdf=dgenbeta(x,a=0.5,b=0.5,c=a,d=b,log=log)
	return(pdf)
}

parcsine=function (x, a=0, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pgenbeta(x,a=0.5,b=0.5,c=a,d=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

vararcsine=function (p, a=0, b=1, log.p=FALSE, lower.tail=TRUE)
{
	var=vargenbeta(p,a=0.5,b=0.5,c=a,d=b,log.p=log.p,lower.tail=lower.tail)
	return(var)
}

esarcsine=function (p, a=0, b=1)
{
	f=function (x) {vararcsine(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Triangular

dtriangular=function (x, a=0, b=2, c=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE&x<a]=0
        pdf[log==FALSE&x>=a&x<=c]=2*(x[x>=a&x<=c]-a)/((b-a)*(c-a))
        pdf[log==FALSE&x>c&x<=b]=2*(b-x[x>c&x<=b])/((b-a)*(b-c))
        pdf[log==FALSE&x>b]=1
        pdf[log==TRUE&x<a]=-Inf
        pdf[log==TRUE&x>=a&x<=c]=log(2*(x[x>=a&x<=c]-a)/((b-a)*(c-a)))
        pdf[log==TRUE&x>c&x<=b]=log(2*(b-x[x>c&x<=b])/((b-a)*(b-c)))
        pdf[log==TRUE&x>b]=0
	return(pdf)
}


ptriangular=function (x, a=0, b=2, c=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE&x<a]=0
        cdf[log.p==FALSE&lower.tail==TRUE&x>=a&x<=c]=(x[x>=a&x<=c]-a)**2/((b-a)*(c-a))
        cdf[log.p==FALSE&lower.tail==TRUE&x>c&x<=b]=1-(b-x[x>c&x<=b])**2/((b-a)*(b-c))
        cdf[log.p==FALSE&lower.tail==TRUE&x>b]=1
        cdf[log.p==FALSE&lower.tail==FALSE&x<a]=1
        cdf[log.p==FALSE&lower.tail==FALSE&x>=a&x<=c]=1-(x[x>=a&x<=c]-a)**2/((b-a)*(c-a))
        cdf[log.p==FALSE&lower.tail==FALSE&x>c&x<=b]=(b-x[x>c&x<=b])**2/((b-a)*(b-c))
        cdf[log.p==FALSE&lower.tail==FALSE&x>b]=0
	cdf[log.p==TRUE&lower.tail==TRUE&x<a]=-Inf
        cdf[log.p==TRUE&lower.tail==TRUE&x>=a&x<=c]=log((x[x>=a&x<=c]-a)**2/((b-a)*(c-a)))
        cdf[log.p==TRUE&lower.tail==TRUE&x>c&x<=b]=log(1-(b-x[x>c&x<=b])**2/((b-a)*(b-c)))
        cdf[log.p==TRUE&lower.tail==TRUE&x>b]=0
        cdf[log.p==TRUE&lower.tail==FALSE&x<a]=0
        cdf[log.p==TRUE&lower.tail==FALSE&x>=a&x<=c]=log(1-(x[x>=a&x<=c]-a)**2/((b-a)*(c-a)))
        cdf[log.p==TRUE&lower.tail==FALSE&x>c&x<=b]=log((b-x[x>c&x<=b])**2/((b-a)*(b-c)))
        cdf[log.p==TRUE&lower.tail==FALSE&x>b]=-Inf
	return(cdf)
}

vartriangular=function (p, a=0, b=2, c=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=p
        var[p<(c-a)/(b-a)]=a+sqrt(p[p<(c-a)/(b-a)]*(b-a)*(c-a))
        var[p>=(c-a)/(b-a)]=b-sqrt((1-p[p>=(c-a)/(b-a)])*(b-a)*(b-c))
	return(var)
}

estriangular=function (p, a=0, b=2, c=1)
{
	f=function (x) {vartriangular(x,a=a,b=b,c=c)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Generalized beta II

dgenbeta2=function (x, a=1, b=1, c=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=c*x**(c-1)*dbeta(x**c,shape1=a,shape2=b)
        pdf[log==TRUE]=dbeta(x**c,shape1=a,shape2=b,log=TRUE)+log(c)+(c-1)*log(x)
	return(pdf)
}

pgenbeta2=function (x, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(x**c,shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}


vargenbeta2=function (p, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(qbeta(p,shape1=a,shape2=b))**(1/c)
	return(var)
}

esgenbeta2=function (p, a=1, b=1, c=1)
{
	f=function (x) {vargenbeta2(x,a=a,b=b,c=c)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}


#Inverse beta

dinvbeta=function (x, a=1, b=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(1+x)**(-2)*dbeta(x/(1+x),shape1=a,shape2=b)
        pdf[log==TRUE]=dbeta(x/(1+x),shape1=a,shape2=b,log=TRUE)-2*log(1+x)
	return(pdf)
}

pinvbeta=function (x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(x/(1+x),shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}


varinvbeta=function (p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(qbeta(p,shape1=a,shape2=b))/(1-qbeta(p,shape1=a,shape2=b))
	return(var)
}

esinvbeta=function (p, a=1, b=1)
{
	f=function (x) {varinvbeta(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Generalized inverse beta

dgeninvbeta=function (x, a=1, c=1, d=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=a*x**(a-1)*(1+x**a)**(-2)*dbeta(x**a/(1+x**a),shape1=c,shape2=d)
        pdf[log==TRUE]=log(a)+(a-1)*log(x)-2*log(1+x**a)+dbeta(x**a/(1+x**a),shape1=c,shape2=d,log=TRUE)
	return(pdf)
}

pgeninvbeta=function (x, a=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(x**a/(1+x**a),shape1=c,shape2=d,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}


vargeninvbeta=function (p, a=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=((qbeta(p,shape1=c,shape2=d))/(1-qbeta(p,shape1=c,shape2=d)))**(1/a)
	return(var)
}

esgeninvbeta=function (p, a=1, c=1, d=1)
{
	f=function (x) {vargeninvbeta(x,a=a,c=c,d=d)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Two sided power

dtsp=function (x, a=1, theta=0.5, log=FALSE)
{
	pdf=x
        pdf[log==FALSE&x<=theta]=a*(x[x<=theta]/theta)**(a-1)
        pdf[log==FALSE&x>theta]=a*((1-x[x>theta])/(1-theta))**(a-1)
        pdf[log==TRUE&x<=theta]=log(a)+(a-1)*log(x[x<=theta]/theta)
        pdf[log==TRUE&x>theta]=log(a)+(a-1)*log((1-x[x>theta])/(1-theta))
	return(pdf)
}

ptsp=function (x, a=1, theta=0.5, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE&x<=theta]=theta*(x[x<=theta]/theta)**a
        cdf[log.p==FALSE&lower.tail==TRUE&x>theta]=1-(1-theta)*((1-x[x>theta])/(1-theta))**a
        cdf[log.p==FALSE&lower.tail==FALSE&x<=theta]=1-theta*(x[x<=theta]/theta)**a
        cdf[log.p==FALSE&lower.tail==FALSE&x>theta]=(1-theta)*((1-x[x>theta])/(1-theta))**a
        cdf[log.p==TRUE&lower.tail==TRUE&x<=theta]=(1-a)*log(theta)+a*log(x[x<=theta])
        cdf[log.p==TRUE&lower.tail==TRUE&x>theta]=log(1-(1-theta)*((1-x[x>theta])/(1-theta))**a)
        cdf[log.p==TRUE&lower.tail==FALSE&x<=theta]=log(1-theta*(x[x<=theta]/theta)**a)
        cdf[log.p==TRUE&lower.tail==FALSE&x>theta]=(1-a)*log(1-theta)+a*log(1-x[x>theta])
	return(cdf)
}

vartsp=function (p, a=1, theta=0.5, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=p
        var[p<=theta]=theta*(p[p<=theta]/theta)**(1/a)
        var[p>theta]=1-(1-theta)*((1-p[p>theta])/(1-theta))**(1/a)
	return(var)
}

estsp=function (p, a=1, theta=0.5)
{
	f=function (x) {vartsp(x,a=a,theta=theta)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Kumaraswamy

dkum=function (x, a=1, b=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=a*b*x**(a-1)*(1-x**a)**(b-1)
        pdf[log==TRUE]=log(a*b)+(a-1)*log(x)+(b-1)*log(1-x**a)
	return(pdf)
}

pkum=function (x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-(1-x**a)**b
        cdf[log.p==FALSE&lower.tail==FALSE]=(1-x**a)**b
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1-(1-x**a)**b)
        cdf[log.p==TRUE&lower.tail==FALSE]=b*log(1-x**a)
	return(cdf)
}

varkum=function (p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(1-(1-p)**(1/b))**(1/a)
	return(var)
}

eskum=function (p, a=1, b=1)
{
	f=function (x) {varkum(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Normal

dnormal=function (x, mu=0, sigma=1, log=FALSE)
{
	pdf=dnorm(x,mean=mu,sd=sigma,log=log)
	return(pdf)
}

pnormal=function (x, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pnorm(x,mean=mu,sd=sigma,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varnormal=function (p, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	var=qnorm(p,mean=mu,sd=sigma,log.p=log.p,lower.tail=lower.tail)
	return(var)
}

esnormal=function (p, mu=0, sigma=1)
{
	f=function (x) {varnormal(x,mu=mu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Kumaraswamy normal

dkumnormal=function (x, mu=0, sigma=1, a=1, b=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=a*b*dnorm(x,mean=mu,sd=sigma)*(pnorm(x,mean=mu,sd=sigma))**(a-1)*(1-(pnorm(x,mean=mu,sd=sigma))**a)**(b-1)
        pdf[log==TRUE]=log(a*b)+dnorm(x,mean=mu,sd=sigma,log=TRUE)+(a-1)*pnorm(x,mean=mu,sd=sigma,log.p=TRUE)+(b-1)*log(1-(pnorm(x,mean=mu,sd=sigma))**a)
	return(pdf)
}

pkumnormal=function (x, mu=0, sigma=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-(1-(pnorm(x,mean=mu,sd=sigma))**a)**b
        cdf[log.p==FALSE&lower.tail==FALSE]=(1-(pnorm(x,mean=mu,sd=sigma))**a)**b
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-(1-(pnorm(x,mean=mu,sd=sigma))**a)**b)
        cdf[log.p==TRUE&lower.tail==FALSE]=b*log(1-(pnorm(x,mean=mu,sd=sigma))**a)
	return(cdf)
}


varkumnormal=function (p, mu=0, sigma=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=mu+sigma*qnorm((1-(1-p)**(1/b))**(1/a))
	return(var)
}

eskumnormal=function (p, mu=0, sigma=1, a=1, b=1)
{
	f=function (x) {varkumnormal(x,mu=mu,sigma=sigma,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Exponential power

dexppower=function (x, mu=0, sigma=1, a=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(exp(-(abs(x-mu))**a*sigma**(-a)/a))/(2*a**(1/a)*sigma*gamma(1+1/a))
        pdf[log==TRUE]=-(abs(x-mu))**a*sigma**(-a)/a-log(2)-(1/a)*log(a)-log(sigma)-lgamma(1+1/a)
	return(pdf)
}

pexppower=function (x, mu=0, sigma=1, a=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE&x<=mu]=0.5-0.5*pgamma((mu-x[x<=mu])**a/(a*sigma**a),shape=1/a)
        cdf[log.p==FALSE&lower.tail==TRUE&x>mu]=0.5+0.5*pgamma((x[x>mu]-mu)**a/(a*sigma**a),shape=1/a)
        cdf[log.p==FALSE&lower.tail==FALSE&x<=mu]=0.5+0.5*pgamma((mu-x[x<=mu])**a/(a*sigma**a),shape=1/a)
        cdf[log.p==FALSE&lower.tail==FALSE&x>mu]=0.5-0.5*pgamma((x[x>mu]-mu)**a/(a*sigma**a),shape=1/a)
        cdf[log.p==TRUE&lower.tail==TRUE&x<=mu]=log(0.5-0.5*pgamma((mu-x[x<=mu])**a/(a*sigma**a),shape=1/a))
        cdf[log.p==TRUE&lower.tail==TRUE&x>mu]=log(0.5+0.5*pgamma((x[x>mu]-mu)**a/(a*sigma**a),shape=1/a))
        cdf[log.p==TRUE&lower.tail==FALSE&x<=mu]=log(0.5+0.5*pgamma((mu-x[x<=mu])**a/(a*sigma**a),shape=1/a))
        cdf[log.p==TRUE&lower.tail==FALSE&x>mu]=log(0.5-0.5*pgamma((x[x>mu]-mu)**a/(a*sigma**a),shape=1/a))
	return(cdf)
}

varexppower=function (p, mu=0, sigma=1, a=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=p
        var[p<=0.5]=mu-a**(1/a)*sigma*(qgamma(1-2*p[p<=0.5],shape=1/a))**(1/a)
        var[p>0.5]=mu+a**(1/a)*sigma*(qgamma(2*p[p>0.5]-1,shape=1/a))**(1/a)
	return(var)
}

esexppower=function (p, mu=0, sigma=1, a=1)
{
	f=function (x) {varexppower(x,mu=mu,sigma=sigma,a=a)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Asymmetric exponential power

daep=function (x, q1=1, q2=1, alpha=0.5, log=FALSE)
{
	K1= 1/(2*q1**(1/q1)*gamma(1+1/q1))
	K2= 1/(2*q2**(1/q2)*gamma(1+1/q2))
	alphastar=alpha*K1/(alpha*K1+(1-alpha)*K2)
	pdf=x
        pdf[log==FALSE&x<=0]=alpha*(exp(-(1/q1)*(abs(x[x<=0]/(2*alphastar)))**q1))/(2*alphastar*q1**(1/q1)*gamma(1+1/q1))
        pdf[log==FALSE&x>0]=(1-alpha)*(exp(-(1/q2)*(abs(x[x>0]/(2-2*alphastar)))**q2))/(2*(1-alphastar)*q2**(1/q2)*gamma(1+1/q2))
        pdf[log==TRUE&x<=0]=log(alpha)-(1/q1)*(abs(x[x<=0]/(2*alphastar)))**q1-log(2)-log(alphastar)-(1/q1)*log(q1)-lgamma(1+1/q1)
        pdf[log==TRUE&x>0]=log(1-alpha)-(1/q2)*(abs(x[x>0]/(2-2*alphastar)))**q2-log(2)-log(1-alphastar)-(1/q2)*log(q2)-lgamma(1+1/q2)
	return(pdf)
}

paep=function (x, q1=1, q2=1, alpha=0.5, log.p=FALSE, lower.tail=TRUE)
{
	K1= 1/(2*q1**(1/q1)*gamma(1+1/q1))
	K2= 1/(2*q2**(1/q2)*gamma(1+1/q2))
	alphastar=alpha*K1/(alpha*K1+(1-alpha)*K2)
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE&x<=0]=alpha-alpha*pgamma((1/q1)*(-x[x<=0]/(2*alphastar))**q1,shape=1/q1)
        cdf[log.p==FALSE&lower.tail==TRUE&x>0]=alpha+(1-alpha)*pgamma((1/q2)*(-x[x>0]/(2-2*alphastar))**q2,shape=1/q2)
        cdf[log.p==FALSE&lower.tail==FALSE&x<=0]=1-alpha+alpha*pgamma((1/q1)*(-x[x<=0]/(2*alphastar))**q1,shape=1/q1)
        cdf[log.p==FALSE&lower.tail==FALSE&x>0]=1-alpha-(1-alpha)*pgamma((1/q2)*(-x[x>0]/(2-2*alphastar))**q2,shape=1/q2)
        cdf[log.p==TRUE&lower.tail==TRUE&x<=0]=log(alpha-alpha*pgamma((1/q1)*(-x[x<=0]/(2*alphastar))**q1,shape=1/q1))
        cdf[log.p==TRUE&lower.tail==TRUE&x>0]=log(alpha+(1-alpha)*pgamma((1/q2)*(-x[x>0]/(2-2*alphastar))**q2,shape=1/q2))
        cdf[log.p==TRUE&lower.tail==FALSE&x<=0]=log(1-alpha+alpha*pgamma((1/q1)*(-x[x<=0]/(2*alphastar))**q1,shape=1/q1))
        cdf[log.p==TRUE&lower.tail==FALSE&x>0]=log(1-alpha-(1-alpha)*pgamma((1/q2)*(-x[x>0]/(2-2*alphastar))**q2,shape=1/q2))
	return(cdf)
}

varaep=function (p, q1=1, q2=1, alpha=0.5, log.p=FALSE, lower.tail=TRUE)
{
	K1= 1/(2*q1**(1/q1)*gamma(1+1/q1))
	K2= 1/(2*q2**(1/q2)*gamma(1+1/q2))
	alphastar=alpha*K1/(alpha*K1+(1-alpha)*K2)
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=p
        var[p<=alpha]=-2*alphastar*(q1*qgamma(1-p[p<=alpha]/alpha,shape=1/q1))**(1/q1)
        var[p>alpha]=-2*(1-alphastar)*(q2*qgamma(1-(1-p[p>alpha])/(1-alpha),shape=1/q2))**(1/q2)
	return(var)
}

esaep=function (p, q1=1, q2=1, alpha=0.5)
{
	f=function (x) {varaep(x,q1=q1,q2=q2,alpha=alpha)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Beta normal


dbetanorm=function (x, mu=0, sigma=1, a=1, b=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=dnorm(x,mean=mu,sd=sigma)*dbeta(pnorm(x,mean=mu,sd=sigma),shape1=a,shape2=b)
        pdf[log==TRUE]=dnorm(x,mean=mu,sd=sigma,log=TRUE)+dbeta(pnorm(x,mean=mu,sd=sigma),shape1=a,shape2=b,log=TRUE)
	return(pdf)
}

pbetanorm=function (x, mu=0, sigma=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(pnorm(x,mean=mu,sd=sigma),shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varbetanorm=function (p, mu=0, sigma=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=mu+sigma*qnorm(qbeta(p,shape1=a,shape2=b))
	return(var)
}

esbetanorm=function (p, mu=0, sigma=1, a=1, b=1)
{
	f=function (x) {varbetanorm(x,mu=mu,sigma=sigma,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}


#Half normal


dhalfnorm=function (x, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=2*dnorm(x,sd=sigma)
        pdf[log==TRUE]=log(2)+dnorm(x,sd=sigma,log=TRUE)
	return(pdf)
}

phalfnorm=function (x, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=2*pnorm(x,sd=sigma)-1
        cdf[log.p==FALSE&lower.tail==FALSE]=2-2*pnorm(x,sd=sigma)
        cdf[log.p==TRUE&lower.tail==TRUE]=log(2*pnorm(x,sd=sigma)-1)
        cdf[log.p==TRUE&lower.tail==FALSE]=log(2-2*pnorm(x,sd=sigma))
	return(cdf)
}

varhalfnorm=function (p, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=sigma*qnorm(0.5*(1+p))
	return(var)
}

eshalfnorm=function (p, sigma=1)
{
	f=function (x) {varhalfnorm(x,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Kumaraswamy half normal

dkumhalfnorm=function (x, sigma=1, a=1, b=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=2*a*b*dnorm(x,sd=sigma)*(2*pnorm(x,sd=sigma)-1)**(a-1)*(1-(2*pnorm(x,sd=sigma)-1)**a)**(b-1)
        pdf[log==TRUE]=log(2*a*b)+dnorm(x,sd=sigma,log=TRUE)+(a-1)*log(2*pnorm(x,sd=sigma)-1)+(b-1)*log(1-(2*pnorm(x,sd=sigma)-1)**a)
	return(pdf)
}

pkumhalfnorm=function (x, sigma=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-(1-(2*pnorm(x,sd=sigma)-1)**a)**b
        cdf[log.p==FALSE&lower.tail==FALSE]=(1-(2*pnorm(x,sd=sigma)-1)**a)**b
        cdf[log.p==TRUE&lower.tail==TRUE]=log(1-(1-(2*pnorm(x,sd=sigma)-1)**a)**b)
        cdf[log.p==TRUE&lower.tail==FALSE]=b*log(1-(2*pnorm(x,sd=sigma)-1)**a)
	return(cdf)
}

varkumhalfnorm=function (p, sigma=1, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=sigma*qnorm(0.5+0.5*(1-(1-p)**(1/b))**(1/a))
	return(var)
}

eskumhalfnorm=function (p, sigma=1, a=1, b=1)
{
	f=function (x) {varkumhalfnorm(x,sigma=sigma,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Student's t

dT=function (x, n=1, log=FALSE)
{
	pdf=dt(x,df=n,log=log)
	return(pdf)
}

pT=function (x, n=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pt(x,df=n,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varT=function (p, n=1, log.p=FALSE, lower.tail=TRUE)
{
	var=qt(p,df=n,log.p=log.p,lower.tail=lower.tail)
	return(var)
}

esT=function (p, n=1)
{
	f=function (x) {varT(x,n=n)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Generalized asymmetric Student's t


dast=function (x, nu1=1, nu2=1, alpha=0.5, log=FALSE)
{
	K1= gamma((nu1+1)/2)/(sqrt(nu1/2)*gamma(nu1/2))
	K2= gamma((nu2+1)/2)/(sqrt(nu2/2)*gamma(nu2/2))
	alphastar=alpha*K1/(alpha*K1+(1-alpha)*K2)
	pdf=x
        pdf[log==FALSE&x<=0]=(alpha/alphastar)*dt(x[x<=0]/(2*alphastar),df=nu1)
        pdf[log==FALSE&x>0]=((1-alpha)/(1-alphastar))*dt(x[x>0]/(2*(1-alphastar)),df=nu2)
        pdf[log==TRUE&x<=0]=log(alpha/alphastar)+dt(x[x<=0]/(2*alphastar),df=nu1,log=TRUE)
        pdf[log==TRUE&x>0]=log((1-alpha)/(1-alphastar))+dt(x[x>0]/(2*(1-alphastar)),df=nu2,log=TRUE)
	return(pdf)
}

past=function (x, nu1=1, nu2=1, alpha=0.5, log.p=FALSE, lower.tail=TRUE)
{
	K1= gamma((nu1+1)/2)/(sqrt(nu1/2)*gamma(nu1/2))
	K2= gamma((nu2+1)/2)/(sqrt(nu2/2)*gamma(nu2/2))
	alphastar=alpha*K1/(alpha*K1+(1-alpha)*K2)
	cdf=x
        mm=x
        MM=x
        for (i in 1:length(x))
        {mm[i]=min(x[i],0)
         MM[i]=max(x[i],0)} 
        cdf[log.p==FALSE&lower.tail==TRUE]=2*alpha*pt(mm/(2*alphastar),df=nu1)-1+alpha+2*(1-alpha)*pt(MM/(2-2*alphastar),df=nu2)
        cdf[log.p==FALSE&lower.tail==FALSE]=1-2*alpha*pt(mm/(2*alphastar),df=nu1)+1-alpha-2*(1-alpha)*pt(MM/(2-2*alphastar),df=nu2)
        cdf[log.p==TRUE&lower.tail==TRUE]=log(2*alpha*pt(mm/(2*alphastar),df=nu1)-1+alpha+2*(1-alpha)*pt(MM/(2-2*alphastar),df=nu2))
        cdf[log.p==TRUE&lower.tail==FALSE]=log(1-2*alpha*pt(mm/(2*alphastar),df=nu1)+1-alpha-2*(1-alpha)*pt(MM/(2-2*alphastar),df=nu2))
	return(cdf)
}


varast=function (p, nu1=1, nu2=1, alpha=0.5, log.p=FALSE, lower.tail=TRUE)
{
	K1= gamma((nu1+1)/2)/(sqrt(nu1/2)*gamma(nu1/2))
	K2= gamma((nu2+1)/2)/(sqrt(nu2/2)*gamma(nu2/2))
	alphastar=alpha*K1/(alpha*K1+(1-alpha)*K2)
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=p
        mm=p
        MM=p
        for (i in 1:length(p))
        {mm[i]=min(p[i],alpha)
         MM[i]=max(p[i],alpha)}
        var=2*alphastar*qt(mm/(2*alpha),df=nu1)+2*(1-alphastar)*qt((MM+1-2*alpha)/(2-2*alpha),df=nu2)
	return(var)
}

esast=function (p, nu1=1, nu2=1, alpha=0.5)
{
	f=function (x) {varast(x,nu1=nu1,nu2=nu2,alpha=alpha)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Half Student's t

dhalfT=function (x, n=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=2*dt(x,df=n)
        pdf[log==TRUE]=log(2)+dt(x,df=n,log=TRUE)
	return(pdf)
}

phalfT=function (x, n=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(x*x/(n+x*x),shape1=0.5,shape2=n/2,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varhalfT=function (p, n=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=sqrt(n*qbeta(p,shape1=0.5,shape2=n/2)/(1-qbeta(p,shape1=0.5,shape2=n/2)))
	return(var)
}

eshalfT=function (p, n=1)
{
	f=function (x) {varhalfT(x,n=n)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Cauchy

dCauchy=function (x, mu=0, sigma=1, log=FALSE)
{
	pdf=dcauchy(x,location=mu,scale=sigma,log=log)
	return(pdf)
}

pCauchy=function (x, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pcauchy(x,location=mu,scale=sigma,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varCauchy=function (p, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	var=qcauchy(p,location=mu,scale=sigma,log.p=log.p,lower.tail=lower.tail)
	return(var)
}

esCauchy=function (p, mu=0, sigma=1)
{
	f=function (x) {varCauchy(x,mu=mu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Log Cauchy

dlogcauchy=function (x, mu=0, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(1/x)*dcauchy(log(x),location=mu,scale=sigma)
        pdf[log==TRUE]=dcauchy(log(x),location=mu,scale=sigma,log=TRUE)-log(x)
	return(pdf)
}

plogcauchy=function (x, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(1/pi)*atan((log(x)-mu)/sigma)
        cdf[log.p==FALSE&lower.tail==FALSE]=1-(1/pi)*atan((log(x)-mu)/sigma)
        cdf[log.p==TRUE&lower.tail==TRUE]=log((1/pi)*atan((log(x)-mu)/sigma))
        cdf[log.p==TRUE&lower.tail==FALSE]=log(1-(1/pi)*atan((log(x)-mu)/sigma))
	return(cdf)
}

varlogcauchy=function (p, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=exp(mu+sigma*tan(pi*p))
	return(var)
}

eslogcauchy=function (p, mu=0, sigma=1)
{
	f=function (x) {varlogcauchy(x,mu=mu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Half Cauchy

dhalfcauchy=function (x, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=2*dcauchy(x,scale=sigma)
        pdf[log==TRUE]=log(2)+dcauchy(x,scale=sigma,log=TRUE)
	return(pdf)
}

phalfcauchy=function (x, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(2/pi)*atan(x/sigma)
        cdf[log.p==FALSE&lower.tail==FALSE]=1-(2/pi)*atan(x/sigma)
        cdf[log.p==TRUE&lower.tail==TRUE]=log((2/pi)*atan(x/sigma))
        cdf[log.p==TRUE&lower.tail==FALSE]=log(1-(2/pi)*atan(x/sigma))
	return(cdf)
}

varhalfcauchy=function (p, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=sigma*tan(pi*p/2)
	return(var)
}

eshalfcauchy=function (p, sigma=1)
{
	f=function (x) {varhalfcauchy(x,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Laplace

dlaplace=function (x, mu=0, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(1/(2*sigma))*exp(-abs(x-mu)/sigma)
        pdf[log==TRUE]=-abs(x-mu)/sigma-log(2*sigma)
	return(pdf)
}

plaplace=function (x, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE&x<mu]=0.5*exp((x[x<mu]-mu)/sigma)
        cdf[log.p==FALSE&lower.tail==TRUE&x>=mu]=1-0.5*exp((mu-x[x>=mu])/sigma)
	cdf[log.p==FALSE&lower.tail==FALSE&x<mu]=1-0.5*exp((x[x<mu]-mu)/sigma)
        cdf[log.p==FALSE&lower.tail==FALSE&x>=mu]=0.5*exp((mu-x[x>=mu])/sigma)
        cdf[log.p==FALSE&lower.tail==TRUE&x<mu]=(x[x<mu]-mu)/sigma-log(2)
        cdf[log.p==FALSE&lower.tail==TRUE&x>=mu]=log(1-0.5*exp((mu-x[x>=mu])/sigma))
	cdf[log.p==FALSE&lower.tail==FALSE&x<mu]=log(1-0.5*exp((x[x<mu]-mu)/sigma))
        cdf[log.p==FALSE&lower.tail==FALSE&x>=mu]=(mu-x[x>=mu])/sigma-log(2)
	return(cdf)
}

varlaplace=function (p, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
        var=p
	var[p<0.5]=mu+sigma*log(2*p[p<0.5])
        var[p>=0.5]=mu-sigma*log(2*(1-p[p>=0.5]))
	return(var)
}

eslaplace=function (p, mu=0, sigma=1)
{
	f=function (x) {varlaplace(x,mu=mu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Poiraud Casanova Thomas Agnan Laplace


dPCTAlaplace=function (x, a=0.5, theta=0, log=FALSE)
{
	pdf=x
        pdf[log==FALSE&x<=theta]=a*(1-a)*exp((1-a)*(x[x<=theta]-theta))
        pdf[log==FALSE&x>theta]=a*(1-a)*exp(a*(theta-x[x>theta]))
        pdf[log==TRUE&x<=theta]=log(a*(1-a))+(1-a)*(x[x<=theta]-theta)
        pdf[log==TRUE&x>theta]=log(a*(1-a))+a*(theta-x[x>theta])
	return(pdf)
}

pPCTAlaplace=function (x, a=0.5, theta=0, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE&x<=theta]=a*exp((1-a)*(x[x<=theta]-theta))
        cdf[log.p==FALSE&lower.tail==TRUE&x>theta]=1-(1-a)*exp(a*(theta-x[x>theta]))
	cdf[log.p==FALSE&lower.tail==FALSE&x<=theta]=1-a*exp((1-a)*(x[x<=theta]-theta))
        cdf[log.p==FALSE&lower.tail==FALSE&x>theta]=(1-a)*exp(a*(theta-x[x>theta]))
	cdf[log.p==TRUE&lower.tail==TRUE&x<=theta]=log(a)+(1-a)*(x[x<=theta]-theta)
        cdf[log.p==TRUE&lower.tail==TRUE&x>theta]=log(1-(1-a)*exp(a*(theta-x[x>theta])))
	cdf[log.p==TRUE&lower.tail==FALSE&x<=theta]=log(1-a*exp((1-a)*(x[x<=theta]-theta)))
        cdf[log.p==TRUE&lower.tail==FALSE&x>theta]=log(1-a)+a*(theta-x[x>theta])
	return(cdf)
}

varPCTAlaplace=function (p, a=0.5, theta=0, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
        var=p
	var[p<=a]=theta+(log(p[p<=a]/a))/(1-a)
        var[p>a]=theta-(log((1-p[p>a])/(1-a)))/a
	return(var)
}

esPCTAlaplace=function (p, a=0.5, theta=0)
{
	f=function (x) {varPCTAlaplace(x,a=a,theta=theta)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Holla Bhattacharya Laplace


dHBlaplace=function (x, a=0.5, theta=0, phi=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE&x<=theta]=a*phi*exp(phi*(x[x<=theta]-theta))
        pdf[log==FALSE&x>theta]=(1-a)*phi*exp(phi*(theta-x[x>theta]))
        pdf[log==TRUE&x<=theta]=log(a*phi)+phi*(x[x<=theta]-theta)
        pdf[log==TRUE&x>theta]=log((1-a)*phi)+phi*(theta-x[x>theta])
	return(pdf)
}

pHBlaplace=function (x, a=0.5, theta=0, phi=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE&x<=theta]=a*exp(phi*(x[x<=theta]-theta))
        cdf[log.p==FALSE&lower.tail==TRUE&x>theta]=1-(1-a)*exp(phi*(theta-x[x>theta]))
	cdf[log.p==FALSE&lower.tail==FALSE&x<=theta]=1-a*exp(phi*(x[x<=theta]-theta))
        cdf[log.p==FALSE&lower.tail==FALSE&x>theta]=(1-a)*exp(phi*(theta-x[x>theta]))
        cdf[log.p==TRUE&lower.tail==TRUE&x<=theta]=log(a)+phi*(x[x<=theta]-theta)
        cdf[log.p==TRUE&lower.tail==TRUE&x>theta]=log(1-(1-a)*exp(phi*(theta-x[x>theta])))
	cdf[log.p==TRUE&lower.tail==FALSE&x<=theta]=log(1-a*exp(phi*(x[x<=theta]-theta)))
        cdf[log.p==TRUE&lower.tail==FALSE&x>theta]=log(1-a)+phi*(theta-x[x>theta])
	return(cdf)
}

varHBlaplace=function (p, a=0.5, theta=0, phi=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
        var=p
	var[p<=a]=theta+(1/phi)*log(p[p<=a]/a)
        var[p>a]=theta-(1/phi)*log((1-p[p>a])/(1-a))
	return(var)
}

esHBlaplace=function (p, a=0.5, theta=0, phi=1)
{
	f=function (x) {varHBlaplace(x,a=a,theta=theta,phi=phi)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#McGill Laplace

dMlaplace=function (x, theta=0, phi=1, psi=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE&x<=theta]=(1/(2*psi))*exp((x[x<=theta]-theta)/psi)
        pdf[log==FALSE&x>theta]=(1/(2*phi))*exp((theta-x[x>theta])/phi)
        pdf[log==TRUE&x<=theta]=(x[x<=theta]-theta)/psi-log(2*psi)
        pdf[log==TRUE&x>theta]=(theta-x[x>theta])/phi-log(2*phi)
	return(pdf)
}

pMlaplace=function (x, theta=0, phi=1, psi=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE&x<=theta]=0.5*exp((x[x<=theta]-theta)/psi)
        cdf[log.p==FALSE&lower.tail==TRUE&x>theta]=1-0.5*exp((theta-x[x>theta])/phi)
	cdf[log.p==FALSE&lower.tail==FALSE&x<=theta]=1-0.5*exp((x[x<=theta]-theta)/psi)
        cdf[log.p==FALSE&lower.tail==FALSE&x>theta]=0.5*exp((theta-x[x>theta])/phi)
	cdf[log.p==TRUE&lower.tail==TRUE&x<=theta]=(x[x<=theta]-theta)/psi-log(2)
        cdf[log.p==TRUE&lower.tail==TRUE&x>theta]=log(1-0.5*exp((theta-x[x>theta])/phi))
	cdf[log.p==TRUE&lower.tail==FALSE&x<=theta]=log(1-0.5*exp((x[x<=theta]-theta)/psi))
        cdf[log.p==TRUE&lower.tail==FALSE&x>theta]=(theta-x[x>theta])/phi-log(2)
	return(cdf)
}

varMlaplace=function (p, theta=0, phi=1, psi=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
        var=p
	var[p<=0.5]=theta+psi*log(2*p[p<=0.5])
        var[p>0.5]=theta-phi*log(2*(1-p[p>0.5]))
	return(var)
}

esMlaplace=function (p, theta=0, phi=1, psi=1)
{
	f=function (x) {varMlaplace(x,theta=theta,phi=phi,psi=psi)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Log Laplace

dloglaplace=function (x, a=1, b=1, delta=0, log=FALSE)
{
	pdf=x
        pdf[log==FALSE&x<=delta]=a*b*x[x<=delta]**(b-1)/(delta**b*(a+b))
        pdf[log==FALSE&x>delta]=a*b*delta**a/(x[x>delta]**(a+1)*(a+b))
        pdf[log==TRUE&x<=delta]=log(a*b)+(b-1)*log(x[x<=delta])-b*log(delta)-log(a+b)
        pdf[log==TRUE&x>delta]=log(a*b)-(a+1)*log(x[x>delta])+a*log(delta)-log(a+b)
	return(pdf)
}

ploglaplace=function (x, a=1, b=1, delta=0, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE&x<=delta]=(a/(a+b))*(x[x<=delta]/delta)**b
        cdf[log.p==FALSE&lower.tail==TRUE&x>delta]=1-(b/(a+b))*(delta/x[x>delta])**a
	cdf[log.p==FALSE&lower.tail==FALSE&x<=delta]=1-(a/(a+b))*(x[x<=delta]/delta)**b
        cdf[log.p==FALSE&lower.tail==FALSE&x>delta]=(b/(a+b))*(delta/x[x>delta])**a
	cdf[log.p==TRUE&lower.tail==TRUE&x<=delta]=log(a/(a+b))+b*log(x[x<=delta]/delta)
        cdf[log.p==TRUE&lower.tail==TRUE&x>delta]=log(1-(b/(a+b))*(delta/x[x>delta])**a)
	cdf[log.p==TRUE&lower.tail==FALSE&x<=delta]=log(1-(a/(a+b))*(x[x<=delta]/delta)**b)
        cdf[log.p==TRUE&lower.tail==FALSE&x>delta]=log(b/(a+b))+a*log(delta/x[x>delta])
	return(cdf)
}

varloglaplace=function (p, a=1, b=1, delta=0, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
        var=p
	var[p<=a/(a+b)]=delta*(p[p<=a/(a+b)]*(a+b)/a)**(1/b)
        var[p>a/(a+b)]=delta*((1-p[p>a/(a+b)])*(a+b)/a)**(-1/a)
	return(var)
}

esloglaplace=function (p, a=1, b=1, delta=0)
{
	f=function (x) {varloglaplace(x,a=a,b=b,delta=delta)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Asymmetric Laplace

dasylaplace=function (x, tau=1, kappa=1, theta=0, log=FALSE)
{
	pdf=x
        pdf[log==FALSE&x>=theta]=sqrt(2)*kappa*(exp(-sqrt(2)*kappa*abs(x[x>=theta]-theta)/tau))/(tau*(1+kappa**2))
        pdf[log==FALSE&x<theta]=sqrt(2)*kappa*(exp(-sqrt(2)*abs(x[x<theta]-theta)/(kappa*tau)))/(tau*(1+kappa**2))
        pdf[log==TRUE&x>=theta]=0.5*log(2)+log(kappa)-sqrt(2)*kappa*abs(x[x>=theta]-theta)/tau-log(tau)-log(1+kappa**2)
        pdf[log==TRUE&x<theta]=0.5*log(2)+log(kappa)-sqrt(2)*abs(x[x<theta]-theta)/(kappa*tau)-log(tau)-log(1+kappa**2)
	return(pdf)
}

pasylaplace=function (x, tau=1, kappa=1, theta=0, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE&x>=theta]=1-(exp(sqrt(2)*kappa*(theta-x[x>=theta])/tau))/(1+kappa**2)
        cdf[log.p==FALSE&lower.tail==TRUE&x<theta]=(kappa**2*exp(sqrt(2)*(x[x<theta]-theta)/(kappa*tau)))/(1+kappa**2)
	cdf[log.p==FALSE&lower.tail==FALSE&x>=theta]=(exp(sqrt(2)*kappa*(theta-x[x>=theta])/tau))/(1+kappa**2)
        cdf[log.p==FALSE&lower.tail==FALSE&x<theta]=1-(kappa**2*exp(sqrt(2)*(x[x<theta]-theta)/(kappa*tau)))/(1+kappa**2)
	cdf[log.p==TRUE&lower.tail==TRUE&x>=theta]=log(1-(exp(sqrt(2)*kappa*(theta-x[x>=theta])/tau))/(1+kappa**2))
        cdf[log.p==TRUE&lower.tail==TRUE&x<theta]=2*log(kappa)+sqrt(2)*(x[x<theta]-theta)/(kappa*tau)-log(1+kappa**2)
	cdf[log.p==TRUE&lower.tail==FALSE&x>=theta]=sqrt(2)*kappa*(theta-x[x>=theta])/tau-log(1+kappa**2)
        cdf[log.p==TRUE&lower.tail==FALSE&x<theta]=log(1-(kappa**2*exp(sqrt(2)*(x[x<theta]-theta)/(kappa*tau)))/(1+kappa**2))
	return(cdf)
}

varasylaplace=function (p, tau=1, kappa=1, theta=0, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
        var=p
	var[p>=kappa**2/(1+kappa**2)]=theta-tau*(log((1-p[p>=kappa**2/(1+kappa**2)])*(1+kappa**2)))/(sqrt(2)*kappa)
        var[p<kappa**2/(1+kappa**2)]=theta+tau*kappa*(log(p[p<kappa**2/(1+kappa**2)]*(1+kappa**(-2))))/(sqrt(2))
	return(var)
}

esasylaplace=function (p, tau=1, kappa=1, theta=0)
{
	f=function (x) {varasylaplace(x,tau=tau,kappa=kappa,theta=theta)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Asymmetric power


dasypower=function (x, a=0.5, lambda=1, delta=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE&x<=0]=delta**(1/lambda)*(exp(-delta*(-x[x<=0]/a)**lambda))/gamma(1+1/lambda)
        pdf[log==FALSE&x>0]=delta**(1/lambda)*(exp(-delta*(x[x>0]/(1-a))**lambda))/gamma(1+1/lambda)
        pdf[log==TRUE&x<=0]=(1/lambda)*log(delta)-delta*(-x[x<=0]/a)**lambda-lgamma(1+1/lambda)
        pdf[log==TRUE&x>0]=(1/lambda)*log(delta)-delta*(x[x>0]/(1-a))**lambda-lgamma(1+1/lambda)
	return(pdf)
}

pasypower=function (x, a=0.5, lambda=1, delta=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE&x<=0]=a-a*dgamma(delta*(-x[x<=0]/a)**lambda,shape=1/lambda)
        cdf[log.p==FALSE&lower.tail==TRUE&x>0]=a-(1-a)*dgamma(delta*(x[x>0]/(1-a))**lambda,shape=1/lambda)
	cdf[log.p==FALSE&lower.tail==FALSE&x<=0]=1-a+a*dgamma(delta*(-x[x<=0]/a)**lambda,shape=1/lambda)
        cdf[log.p==FALSE&lower.tail==FALSE&x>0]=1-a+(1-a)*dgamma(delta*(x[x>0]/(1-a))**lambda,shape=1/lambda)
	cdf[log.p==TRUE&lower.tail==TRUE&x<=0]=log(a-a*dgamma(delta*(-x[x<=0]/a)**lambda,shape=1/lambda))
        cdf[log.p==TRUE&lower.tail==TRUE&x>0]=log(a-(1-a)*dgamma(delta*(x[x>0]/(1-a))**lambda,shape=1/lambda))
	cdf[log.p==TRUE&lower.tail==FALSE&x<=0]=log(1-a+a*dgamma(delta*(-x[x<=0]/a)**lambda,shape=1/lambda))
        cdf[log.p==TRUE&lower.tail==FALSE&x>0]=log(1-a+(1-a)*dgamma(delta*(x[x>0]/(1-a))**lambda,shape=1/lambda))
	return(cdf)
}

varasypower=function (p, a=0.5, lambda=1, delta=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
        var=p
	var[p<=a]=-a*((qgamma(1-p[p<=a]/a,shape=1/lambda))**(1/lambda))/delta**(1/lambda)
        var[p>a]=(1-a)*((qgamma(1-(1-p[p>a])/(1-a),shape=1/lambda))**(1/lambda))/delta**(1/lambda)
	return(var)
}

esasypower=function (p, a=0.5, lambda=1, delta=1)
{
	f=function (x) {varasypower(x,a=a,lambda=lambda,delta=delta)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Logistic

dlogistic=function (x, mu=0, sigma=1, log=FALSE)
{
	pdf=dlogis(x,location=mu,scale=sigma,log=log)
	return(pdf)
}

plogistic=function (x, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=plogis(x,location=mu,scale=sigma,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varlogistic=function (p, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	var=qlogis(p,location=mu,scale=sigma,log.p=log.p,lower.tail=lower.tail)
	return(var)
}

eslogistic=function (p, mu=0, sigma=1)
{
	f=function (x) {varlogistic(x,mu=mu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Hyperbolic secant

dsecant=function (x, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=0.5/cosh(pi*x/2)
        pdf[log==TRUE]=-log(2)-log(cosh(pi*x/2))
	return(pdf)
}

psecant=function (x, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(2/pi)*atan(exp(pi*x/2))
        cdf[log.p==FALSE&lower.tail==FALSE]=1-(2/pi)*atan(exp(pi*x/2))
	cdf[log.p==TRUE&lower.tail==TRUE]=log(2)-log(pi)+log(atan(exp(pi*x/2)))
	cdf[log.p==TRUE&lower.tail==FALSE]=log(1-(2/pi)*atan(exp(pi*x/2)))
	return(cdf)
}

varsecant=function (p, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(2/pi)*log(tan(pi*p/2))
	return(var)
}

essecant=function (p)
{
	f=function (x) {varsecant(x)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Generalized logistic


dgenlogis=function (x, a=1, mu=0, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(a/sigma)*exp(-(x-mu)/sigma)*(1+exp(-(x-mu)/sigma))**(-1-a)
        pdf[log==TRUE]=log(a)-log(sigma)-(x-mu)/sigma-(1+a)*log(1+exp(-(x-mu)/sigma))
	return(pdf)
}

pgenlogis=function (x, a=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(1+exp(-(x-mu)/sigma))**(-a)
        cdf[log.p==FALSE&lower.tail==FALSE]=1-(1+exp(-(x-mu)/sigma))**(-a)
	cdf[log.p==TRUE&lower.tail==TRUE]=-a*log(1+exp(-(x-mu)/sigma))
	cdf[log.p==TRUE&lower.tail==FALSE]=log(1-(1+exp(-(x-mu)/sigma))**(-a))
	return(cdf)
}

vargenlogis=function (p, a=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=mu-sigma*log(p**(-1/a)-1)
	return(var)
}

esgenlogis=function (p, a=1, mu=0, sigma=1)
{
	f=function (x) {vargenlogis(x,a=a,mu=mu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Generalized logistic III

dgenlogis3=function (x, alpha=1, mu=0, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=exp(alpha*(x-mu)/sigma)*((1+exp((x-mu)/sigma))**(-1-alpha))/(sigma*beta(alpha,alpha))
        pdf[log==TRUE]=alpha*(x-mu)/sigma-(1+alpha)*log(1+exp((x-mu)/sigma))-log(sigma)-lbeta(alpha,alpha)
	return(pdf)
}

pgenlogis3=function (x, alpha=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(1/(1+exp((mu-x)/sigma)),shape1=alpha,shape2=alpha,log.p=FALSE,lower.tail=TRUE)
	return(cdf)
}

vargenlogis3=function (p, alpha=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=mu-sigma*log(1/qbeta(p,shape1=alpha,shape2=alpha)-1)
	return(var)
}

esgenlogis3=function (p, alpha=1, mu=0, sigma=1)
{
	f=function (x) {vargenlogis3(x,alpha=alpha,mu=mu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Generalized logistic IV

dgenlogis4=function (x, a=1, alpha=1, mu=0, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=exp(-alpha*(x-mu)/sigma)*((1+exp(-(x-mu)/sigma))**(-a-alpha))/(sigma*beta(a,alpha))
        pdf[log==TRUE]=-alpha*(x-mu)/sigma-(a+alpha)*log(1+exp(-(x-mu)/sigma))-log(sigma)-lbeta(a,alpha)
	return(pdf)
}

pgenlogis4=function (x, a=1, alpha=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(1/(1+exp((mu-x)/sigma)),shape1=alpha,shape2=a,log.p=FALSE,lower.tail=TRUE)
	return(cdf)
}

vargenlogis4=function (p, a=1, alpha=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=mu-sigma*log(1/qbeta(p,shape1=alpha,shape2=a)-1)
	return(var)
}

esgenlogis4=function (p, a=1, alpha=1, mu=0, sigma=1)
{
	f=function (x) {vargenlogis4(x,a=a,alpha=alpha,mu=mu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Half logistic

dhalflogis=function (x, lambda=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=2*lambda*exp(-lambda*x)*(1+exp(-lambda*x))**(-2)
        pdf[log==TRUE]=log(2*lambda)-lambda*x-2*log(1+exp(-lambda*x))
	return(pdf)
}

phalflogis=function (x, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(1-exp(-lambda*x))/(1+exp(-lambda*x))
        cdf[log.p==FALSE&lower.tail==FALSE]=1-(1-exp(-lambda*x))/(1+exp(-lambda*x))
	cdf[log.p==TRUE&lower.tail==TRUE]=log((1-exp(-lambda*x))/(1+exp(-lambda*x)))
	cdf[log.p==TRUE&lower.tail==FALSE]=log(1-(1-exp(-lambda*x))/(1+exp(-lambda*x)))
	return(cdf)
}

varhalflogis=function (p, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=-(1/lambda)*log((1-p)/(1+p))
	return(var)
}

eshalflogis=function (p, lambda=1)
{
	f=function (x) {varhalflogis(x,lambda=lambda)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Log logistic

dloglogis=function (x, a=1, b=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=b*a**b*x**(b-1)*(a**b+x**b)**(-2)
        pdf[log==TRUE]=log(b)+b*log(a)+(b-1)*log(x)-2*log(a**b+x**b)
	return(pdf)
}

ploglogis=function (x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=x**b/(a**b+x**b)
        cdf[log.p==FALSE&lower.tail==FALSE]=a**b/(a**b+x**b)
	cdf[log.p==TRUE&lower.tail==TRUE]=b*log(x)-log(a**b+x**b)
	cdf[log.p==TRUE&lower.tail==FALSE]=b*log(a)-log(a**b+x**b)
	return(cdf)
}

varloglogis=function (p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=a*(p/(1-p))**(1/b)
	return(var)
}

esloglogis=function (p, a=1, b=1)
{
	f=function (x) {varloglogis(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Kumaraswamy log logistic

dkumloglogis=function (x, a=1, b=1, alpha=1, beta=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=a*b*beta*alpha**beta*x**(a*beta-1)*(alpha**beta+x**beta)**(-a-1)*(1-x**(a*beta)*(alpha**beta+x**beta)**(-a))**(b-1)
        pdf[log==TRUE]=log(a*b*beta)+beta*log(alpha)+(a*beta-1)*log(x)-(a+1)*log(alpha**beta+x**beta)+(b-1)*log(1-x**(a*beta)*(alpha**beta+x**beta)**(-a))
	return(pdf)
}

pkumloglogis=function (x, a=1, b=1, alpha=1, beta=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(1-x**(a*beta)*(alpha**beta+x**beta)**(-a))**b
        cdf[log.p==FALSE&lower.tail==FALSE]=1-(1-x**(a*beta)*(alpha**beta+x**beta)**(-a))**b
	cdf[log.p==TRUE&lower.tail==TRUE]=b*log(1-x**(a*beta)*(alpha**beta+x**beta)**(-a))
	cdf[log.p==TRUE&lower.tail==FALSE]=log(1-(1-x**(a*beta)*(alpha**beta+x**beta)**(-a))**b)
	return(cdf)
}

varkumloglogis=function (p, a=1, b=1, alpha=1, beta=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=alpha*((1-(1-p)**(1/b))**(1/a)-1)**(-1/beta)
	return(var)
}

eskumloglogis=function (p, a=1, b=1, alpha=1, beta=1)
{
	f=function (x) {varkumloglogis(x,a=a,b=b,alpha=alpha,beta=beta)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Exponentiated logistic

dexplogis=function (x, a=1, b=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(a/b)*exp(-x/b)*(1+exp(-x/b))**(-a-1)
        pdf[log==TRUE]=log(a)-log(b)-x/b-(a+1)*log(1+exp(-x/b))
	return(pdf)
}

pexplogis=function (x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(1+exp(-x/b))**(-a)
        cdf[log.p==FALSE&lower.tail==FALSE]=1-(1+exp(-x/b))**(-a)
	cdf[log.p==TRUE&lower.tail==TRUE]=-a*log(1+exp(-x/b))
	cdf[log.p==TRUE&lower.tail==FALSE]=log(1-(1+exp(-x/b))**(-a))
	return(cdf)
}

varexplogis=function (p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=-b*log(p**(-1/a)-1)
	return(var)
}

esexplogis=function (p, a=1, b=1)
{
	f=function (x) {varexplogis(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Hosking logistic

dHlogis=function (x, k=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(1-k*x)**(1/k-1)*(1+(1-k*x)**(1/k))**(-2)
        pdf[log==TRUE]=(1/k-1)*log(1-k*x)-2*log(1+(1-k*x)**(1/k))
	return(pdf)
}

pHlogis=function (x, k=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(1+(1-k*x)**(1/k))**(-1)
        cdf[log.p==FALSE&lower.tail==FALSE]=(1-k*x)**(1/k)*(1+(1-k*x)**(1/k))**(-1)
	cdf[log.p==TRUE&lower.tail==TRUE]=-log(1+(1-k*x)**(1/k))
	cdf[log.p==TRUE&lower.tail==FALSE]=(1/k)*log(1-k*x)-log(1+(1-k*x)**(1/k))
	return(cdf)
}

varHlogis=function (p, k=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(1/k)*(1-(1/p-1)**k)
	return(var)
}

esHlogis=function (p, k=1)
{
	f=function (x) {varHlogis(x,k=k)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Lognormal

dlognorm=function (x, mu=0, sigma=1, log=FALSE)
{
	pdf=dlnorm(x,meanlog=mu,sdlog=sigma,log=log)
	return(pdf)
}

plognorm=function (x, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=plnorm(x,meanlog=mu,sdlog=sigma,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varlognorm=function (p, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	var=qlnorm(p,meanlog=mu,sdlog=sigma,log.p=log.p,lower.tail=lower.tail)
	return(var)
}

eslognorm=function (p, mu=0, sigma=1)
{
	f=function (x) {varlognorm(x,mu=mu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Beta lognormal

dbetalognorm=function (x, a=1, b=1, mu=0, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=dlnorm(x,meanlog=mu,sdlog=sigma)*dbeta(plnorm(x,meanlog=mu,sdlog=sigma),shape1=a,shape2=b)
        pdf[log==TRUE]=dlnorm(x,meanlog=mu,sdlog=sigma,log=TRUE)+dbeta(plnorm(x,meanlog=mu,sdlog=sigma),shape1=a,shape2=b,log=TRUE)
	return(pdf)
}

pbetalognorm=function (x, a=1, b=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(plnorm(x,meanlog=mu,sdlog=sigma),shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varbetalognorm=function (p, a=1, b=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=exp(mu+sigma*qnorm(qbeta(p,shape1=a,shape2=b)))
	return(var)
}

esbetalognorm=function (p, a=1, b=1, mu=0, sigma=1)
{
	f=function (x) {varbetalognorm(x,a=a,b=b,mu=mu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Burr

dburr=function (x, a=1, b=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=b*a**b*x**(-b-1)*(1+(x/a)**(-b))**(-2)
        pdf[log==TRUE]=log(b)+b*log(a)-(b+1)*log(x)-2*log(1+(x/a)**(-b))
	return(pdf)
}

pburr=function (x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(1+(x/a)**(-b))**(-1)
        cdf[log.p==FALSE&lower.tail==FALSE]=(1+(x/a)**b)**(-1)
	cdf[log.p==TRUE&lower.tail==TRUE]=-log(1+(x/a)**(-b))
	cdf[log.p==TRUE&lower.tail==FALSE]=-log(1+(x/a)**b)
	return(cdf)
}

varburr=function (p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=a*p**(1/b)*(1-p)**(-1/b)
	return(var)
}

esburr=function (p, a=1, b=1)
{
	f=function (x) {varburr(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Beta Burr

dbetaburr=function (x, a=1, b=1, c=1, d=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=b*a**b*x**(-b-1)*(1+(x/a)**(-b))**(-2)*dbeta((1+(x/a)**(-b))**(-1),shape1=c,shape2=d)
        pdf[log==TRUE]=log(b)+b*log(a)-(b+1)*log(x)-2*log(1+(x/a)**(-b))+dbeta((1+(x/a)**(-b))**(-1),shape1=c,shape2=d,log=TRUE)
	return(pdf)
}


pbetaburr=function (x, a=1, b=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf=pbeta((1+(x/a)**(-b))**(-1),shape1=c,shape2=d,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varbetaburr=function (p, a=1, b=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=a*(qbeta(p,shape1=c,shape2=d)/(1-qbeta(p,shape1=c,shape2=d)))**(1/b)
	return(var)
}

esbetaburr=function (p, a=1, b=1, c=1, d=1)
{
	f=function (x) {varbetaburr(x,a=a,b=b,c=c,d=d)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Burr XII

dburr7=function (x, k=1, c=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=k*c*x**(c-1)*(1+x**c)**(-k-1)
        pdf[log==TRUE]=log(k*c)+(c-1)*log(x)-(k+1)*log(1+x**c)
	return(pdf)
}

pburr7=function (x, k=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-(1+x**c)**(-k)
        cdf[log.p==FALSE&lower.tail==FALSE]=(1+x**c)**(-k)
	cdf[log.p==TRUE&lower.tail==TRUE]=-log(1-(1+x**c)**(-k))
	cdf[log.p==TRUE&lower.tail==FALSE]=-k*log(1+x**c)
	return(cdf)
}

varburr7=function (p, k=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=((1-p)**(-1/k)-1)**(1/c)
	return(var)
}

esburr7=function (p, k=1, c=1)
{
	f=function (x) {varburr7(x,k=k,c=c)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Kumaraswamy Burr XII

dkumburr7=function (x, a=1, b=1, k=1, c=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=a*b*k*c*x**(c-1)*(1+x**c)**(-k-1)*(1-(1+x**c)**(-k))**(a-1)*(1-(1-(1+x**c)**(-k))**a)**(b-1)
        pdf[log==TRUE]=log(a*b*k*c)+(c-1)*log(x)-(k+1)*log(1+x**c)+(a-1)*log(1-(1+x**c)**(-k))+(b-1)*log(1-(1-(1+x**c)**(-k))**a)
	return(pdf)
}

pkumburr7=function (x, a=1, b=1, k=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-(1-(1-(1+x**c)**(-k))**a)**b
        cdf[log.p==FALSE&lower.tail==FALSE]=(1-(1-(1+x**c)**(-k))**a)**b
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1-(1-(1-(1+x**c)**(-k))**a)**b)
	cdf[log.p==TRUE&lower.tail==FALSE]=b*log(1-(1-(1+x**c)**(-k))**a)
	return(cdf)
}

varkumburr7=function (p, a=1, b=1, k=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=((1-(1-(1-p)**(1/b))**(1/a))**(-1/k)-1)**(1/c)
	return(var)
}

eskumburr7=function (p, a=1, b=1, k=1, c=1)
{
	f=function (x) {varkumburr7(x,a=a,b=b,k=k,c=c)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Beta Burr XII

dbetaburr7=function (x, a=1, b=1, c=1, k=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=c*k*x**(c-1)*(1+x**c)**(-k-1)*dbeta(1-(1+x**c)**(-k),shape1=a,shape2=b)
        pdf[log==TRUE]=log(c*k)+(c-1)*log(x)-(k+1)*log(1+x**c)+dbeta(1-(1+x**c)**(-k),shape1=a,shape2=b,log=TRUE)
	return(pdf)
}


pbetaburr7=function (x, a=1, b=1, c=1, k=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf=pbeta(1-(1+x**c)**(-k),shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varbetaburr7=function (p, a=1, b=1, c=1, k=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=((1-qbeta(p,shape1=a,shape2=b))**(-1/k)-1)**(1/c)
	return(var)
}

esbetaburr7=function (p, a=1, b=1, c=1, k=1)
{
	f=function (x) {varbetaburr7(x,a=a,b=b,c=c,k=k)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Dagum


ddagum=function (x, a=1, b=1, c=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=a*c*b**a*x**(a*c-1)*(x**a+b**a)**(-c-1)
        pdf[log==TRUE]=log(a*c)+a*log(b)+(a*c-1)*log(x)-(c+1)*log(x**a+b**a)
	return(pdf)
}

pdagum=function (x, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(1+(b/x)**a)**(-c)
        cdf[log.p==FALSE&lower.tail==FALSE]=1-(1+(b/x)**a)**(-c)
	cdf[log.p==TRUE&lower.tail==TRUE]=-c*log(1+(b/x)**a)
	cdf[log.p==TRUE&lower.tail==FALSE]=log(1-(1+(b/x)**a)**(-c))
	return(cdf)
}

vardagum=function (p, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=b*(1-p**(-1/c))**(-1/a)
	return(var)
}

esdagum=function (p, a=1, b=1, c=1)
{
	f=function (x) {vardagum(x,a=a,b=b,c=c)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Lomax

dlomax=function (x, a=1, lambda=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(a/lambda)*(1+x/lambda)**(-a-1)
        pdf[log==TRUE]=log(a/lambda)-(a+1)*log(1+x/lambda)
	return(pdf)
}

plomax=function (x, a=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-(1+x/lambda)**(-a)
        cdf[log.p==FALSE&lower.tail==FALSE]=(1+x/lambda)**(-a)
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1-(1+x/lambda)**(-a))
	cdf[log.p==TRUE&lower.tail==FALSE]=-a*log(1+x/lambda)
	return(cdf)
}

varlomax=function (p, a=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=lambda*((1-p)**(-1/a)-1)
	return(var)
}

eslomax=function (p, a=1, lambda=1)
{
	f=function (x) {varlomax(x,a=a,lambda=lambda)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Beta Lomax

dbetalomax=function (x, a=1, b=1, alpha=1, lambda=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(alpha/lambda)*(1+x/lambda)**(-alpha-1)*dbeta(1-(1+x/lambda)**(-alpha),shape1=a,shape2=b)
        pdf[log==TRUE]=log(alpha/lambda)-(alpha+1)*log(1+x/lambda)+dbeta(1-(1+x/lambda)**(-alpha),shape1=a,shape2=b,log=TRUE)
	return(pdf)
}


pbetalomax=function (x, a=1, b=1, alpha=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(1-(1+x/lambda)**(-alpha),shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varbetalomax=function (p, a=1, b=1, alpha=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=lambda*((1-qbeta(p,shape1=a,shape2=b))**(-1/alpha)-1)
	return(var)
}

esbetalomax=function (p, a=1, b=1, alpha=1, lambda=1)
{
	f=function (x) {varbetalomax(x,a=a,b=b,alpha=alpha,lambda=lambda)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Gumbel

dgumbel=function (x, mu=0, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(1/sigma)*exp(-(x-mu)/sigma)*exp(-exp(-(x-mu)/sigma))
        pdf[log==TRUE]=-log(sigma)-(x-mu)/sigma-exp(-(x-mu)/sigma)
	return(pdf)
}

pgumbel=function (x, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=exp(-exp(-(x-mu)/sigma))
        cdf[log.p==FALSE&lower.tail==FALSE]=1-exp(-exp(-(x-mu)/sigma))
	cdf[log.p==TRUE&lower.tail==TRUE]=-exp(-(x-mu)/sigma)
	cdf[log.p==TRUE&lower.tail==FALSE]=log(1-exp(-exp(-(x-mu)/sigma)))
	return(cdf)
}

vargumbel=function (p, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=mu-sigma*log(-log(p))
	return(var)
}

esgumbel=function (p, mu=0, sigma=1)
{
	f=function (x) {vargumbel(x,mu=mu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Kumaraswamy Gumbel

dkumgumbel=function (x, a=1, b=1, mu=0, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(a*b/sigma)*exp(-(x-mu)/sigma)*exp(-a*exp(-(x-mu)/sigma))*(1-exp(-a*exp(-(x-mu)/sigma)))**(b-1)
        pdf[log==TRUE]=-log(a*b/sigma)-(x-mu)/sigma-a*exp(-(x-mu)/sigma)+(b-1)*log(1-exp(-a*exp(-(x-mu)/sigma)))
	return(pdf)
}

pkumgumbel=function (x, a=1, b=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-(1-exp(-a*exp(-(x-mu)/sigma)))**b
        cdf[log.p==FALSE&lower.tail==FALSE]=(1-exp(-a*exp(-(x-mu)/sigma)))**b
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1-(1-exp(-a*exp(-(x-mu)/sigma)))**b)
	cdf[log.p==TRUE&lower.tail==FALSE]=b*log(1-exp(-a*exp(-(x-mu)/sigma)))
	return(cdf)
}

varkumgumbel=function (p, a=1, b=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=mu-sigma*log(-(1/a)*log(1-(1-p)**(1/b)))
	return(var)
}

eskumgumbel=function (p, a=1, b=1, mu=0, sigma=1)
{
	f=function (x) {varkumgumbel(x,a=a,b=b,mu=mu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Beta Gumbel

dbetagumbel=function (x, a=1, b=1, mu=0, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=(1/sigma)*exp(-(x-mu)/sigma)*exp(-exp(-(x-mu)/sigma))*dbeta(exp(-exp(-(x-mu)/sigma)),shape1=a,shape2=b)
        pdf[log==TRUE]=-log(sigma)-(x-mu)/sigma-exp(-(x-mu)/sigma)+dbeta(exp(-exp(-(x-mu)/sigma)),shape1=a,shape2=b,log=TRUE)
	return(pdf)
}

pbetagumbel=function (x, a=1, b=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(exp(-exp(-(x-mu)/sigma)),shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varbetagumbel=function (p, a=1, b=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=mu-sigma*log(-log(qbeta(p,shape1=a,shape2=b)))
	return(var)
}

esbetagumbel=function (p, a=1, b=1, mu=0, sigma=1)
{
	f=function (x) {varbetagumbel(x,a=a,b=b,mu=mu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Gumbel II

dgumbel2=function (x, a=1, b=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=a*b*x**(-a-1)*exp(-b*x**(-a))
        pdf[log==TRUE]=log(a*b)-(a+1)*log(x)-b*x**(-a)
	return(pdf)
}

pgumbel2=function (x, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-exp(-b*x**(-a))
        cdf[log.p==FALSE&lower.tail==FALSE]=exp(-b*x**(-a))
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1-exp(-b*x**(-a)))
	cdf[log.p==TRUE&lower.tail==FALSE]=-b*x**(-a)
	return(cdf)
}

vargumbel2=function (p, a=1, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=b**(1/a)*(-log(1-p))**(-1/a)
	return(var)
}

esgumbel2=function (p, a=1, b=1)
{
	f=function (x) {vargumbel2(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Beta Gumbel II


dbetagumbel2=function (x, a=1, b=1, c=1, d=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=a*b*x**(-a-1)*exp(-b*x**(-a))*dbeta(1-exp(-b*x**(-a)),shape1=c,shape2=d)
        pdf[log==TRUE]=log(a*b)-(a+1)*log(x)-b*x**(-a)+dbeta(1-exp(-b*x**(-a)),shape1=c,shape2=d,log=TRUE)
	return(pdf)
}


pbetagumbel2=function (x, a=1, b=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(1-exp(-b*x**(-a)),shape1=c,shape2=d,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varbetagumbel2=function (p, a=1, b=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=b**(1/a)*(-log(1-qbeta(p,shape1=c,shape2=d)))**(-1/a)
	return(var)
}

esbetagumbel2=function (p, a=1, b=1, c=1, d=1)
{
	f=function (x) {varbetagumbel2(x,a=a,b=b,c=c,d=d)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Frechet

dfrechet=function (x, alpha=1, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=alpha*sigma**alpha*x**(-alpha-1)*exp(-(sigma/x)**alpha)
        pdf[log==TRUE]=log(alpha)+alpha*log(sigma)-(alpha+1)*log(x)-(sigma/x)**alpha
	return(pdf)
}

pfrechet=function (x, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=exp(-(sigma/x)**alpha)
        cdf[log.p==FALSE&lower.tail==FALSE]=1-exp(-(sigma/x)**alpha)
	cdf[log.p==TRUE&lower.tail==TRUE]=-(sigma/x)**alpha
	cdf[log.p==TRUE&lower.tail==FALSE]=log(1-exp(-(sigma/x)**alpha))
	return(cdf)
}

varfrechet=function (p, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=sigma*(-log(p))**(-1/alpha)
	return(var)
}

esfrechet=function (p, alpha=1, sigma=1)
{
	f=function (x) {varfrechet(x,alpha=alpha,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Beta Frechet

dbetafrechet=function (x, a=1, b=1, alpha=1, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=alpha*sigma**alpha*x**(-alpha-1)*exp(-(sigma/x)**alpha)*dbeta(exp(-(sigma/x)**alpha),shape1=a,shape2=b)
        pdf[log==TRUE]=log(alpha)+alpha*log(sigma)-(alpha+1)*log(x)-(sigma/x)**alpha+dbeta(exp(-(sigma/x)**alpha),shape1=a,shape2=b,log=TRUE)
	return(pdf)
}


pbetafrechet=function (x, a=1, b=1, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(exp(-(sigma/x)**alpha),shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varbetafrechet=function (p, a=1, b=1, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=sigma*(-log(qbeta(p,shape1=a,shape2=b)))**(-1/alpha)
	return(var)
}

esbetafrechet=function (p, a=1, b=1, alpha=1, sigma=1)
{
	f=function (x) {varbetafrechet(x,a=a,b=b,alpha=alpha,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Weibull

dWeibull=function (x, alpha=1, sigma=1, log=FALSE)
{
	pdf=dweibull(x,shape=alpha,scale=sigma,log=log)
	return(pdf)
}

pWeibull=function (x, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pweibull(x,shape=alpha,scale=sigma,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varWeibull=function (p, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	var=qweibull(p,shape=alpha,scale=sigma,log.p=log.p,lower.tail=lower.tail)
	return(var)
}

esWeibull=function (p, alpha=1, sigma=1)
{
	f=function (x) {varWeibull(x,alpha=alpha,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Kumaraswamy Weibull

dkumweibull=function (x, a=1, b=1, alpha=1, sigma=1, log=FALSE)
{
			pdf=x
	pdf[log==FALSE]=a*b*dweibull(x,shape=alpha,scale=sigma)*(1-pweibull(x,shape=alpha,scale=sigma))**(a-1)*(1-(1-pweibull(x,shape=alpha,scale=sigma))**a)**(b-1)
        pdf[log==TRUE]=log(a*b)+dweibull(x,shape=alpha,scale=sigma,log=TRUE)+(a-1)*log(1-pweibull(x,shape=alpha,scale=sigma))+(b-1)*log(1-(1-pweibull(x,shape=alpha,scale=sigma))**a)
	return(pdf)
}

pkumweibull=function (x, a=1, b=1, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-(1-(1-pweibull(x,shape=alpha,scale=sigma))**a)**b
        cdf[log.p==FALSE&lower.tail==FALSE]=(1-(1-pweibull(x,shape=alpha,scale=sigma))**a)**b
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1-(1-(1-pweibull(x,shape=alpha,scale=sigma))**a)**b)
	cdf[log.p==TRUE&lower.tail==FALSE]=b*log(1-(1-pweibull(x,shape=alpha,scale=sigma))**a)
	return(cdf)
}

varkumweibull=function (p, a=1, b=1, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=sigma*(-log(1-(1-(1-p)**(1/b))**(1/a)))**(1/alpha)
	return(var)
}

eskumweibull=function (p, a=1, b=1, alpha=1, sigma=1)
{
	f=function (x) {varkumweibull(x,a=a,b=b,alpha=alpha,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Logistic Rayleigh

dlogisrayleigh=function (x, a=1, lambda=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=a*lambda*x*exp(lambda*x*x/2)*(exp(lambda*x*x/2)-1)**(a-1)*(1+(exp(lambda*x*x/2)-1)**a)**(-2)
        pdf[log==TRUE]=log(a*lambda)+log(x)+lambda*x*x/2+(a-1)*log(exp(lambda*x*x/2)-1)-2*log(1+(exp(lambda*x*x/2)-1)**a)
	return(pdf)
}

plogisrayleigh=function (x, a=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-1/(1+(exp(lambda*x*x/2)-1)**a)
        cdf[log.p==FALSE&lower.tail==FALSE]=1/(1+(exp(lambda*x*x/2)-1)**a)
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1-1/(1+(exp(lambda*x*x/2)-1)**a))
	cdf[log.p==TRUE&lower.tail==FALSE]=-log(1+(exp(lambda*x*x/2)-1)**a)
	return(cdf)
}

varlogisrayleigh=function (p, a=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=sqrt(2/lambda)*sqrt(log(1+(p/(1-p))**(1/a)))
	return(var)
}

eslogisrayleigh=function (p, a=1, lambda=1)
{
	f=function (x) {varlogisrayleigh(x,a=a,lambda=lambda)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Marshall Olkin Weibull

dmoweibull=function (x, a=1, b=1, lambda=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=b*lambda**b*x**(b-1)*exp((lambda*x)**b)*(exp((lambda*x)**b)-1+a)**(-2)
        pdf[log==TRUE]=log(b)+b*log(lambda)+(b-1)*log(x)+(lambda*x)**b-2*log(exp((lambda*x)**b)-1+a)
	return(pdf)
}

pmoweibull=function (x, a=1, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-1/(exp((lambda*x)**b)-1+a)
        cdf[log.p==FALSE&lower.tail==FALSE]=1/(exp((lambda*x)**b)-1+a)
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1-1/(exp((lambda*x)**b)-1+a))
	cdf[log.p==TRUE&lower.tail==FALSE]=-log(exp((lambda*x)**b)-1+a)
	return(cdf)
}

varmoweibull=function (p, a=1, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(1/lambda)*(log(1/(1-p)+1-a))**(1/b)
	return(var)
}

esmoweibull=function (p, a=1, b=1, lambda=1)
{
	f=function (x) {varmoweibull(x,a=a,b=b,lambda=lambda)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Beta Weibull

dbetaweibull=function (x, a=1, b=1, alpha=1, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=dweibull(x,shape=alpha,scale=sigma)*dbeta(pweibull(x,shape=alpha,scale=sigma),shape1=a,shape2=b)
        pdf[log==TRUE]=dweibull(x,shape=alpha,scale=sigma,log=TRUE)+dbeta(pweibull(x,shape=alpha,scale=sigma),shape1=a,shape2=b,log=TRUE)
	return(pdf)
}


pbetaweibull=function (x, a=1, b=1, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pbeta(pweibull(x,shape=alpha,scale=sigma),shape1=a,shape2=b,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varbetaweibull=function (p, a=1, b=1, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=sigma*(-log(1-qbeta(p,shape1=a,shape2=b)))**(1/alpha)
	return(var)
}

esbetaweibull=function (p, a=1, b=1, alpha=1, sigma=1)
{
	f=function (x) {varbetaweibull(x,a=a,b=b,alpha=alpha,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Double Weibull

ddweibull=function (x, c=1, mu=0, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=0.5*dweibull(x-mu,shape=c,scale=sigma)
        pdf[log==TRUE]=dweibull(x-mu,shape=c,scale=sigma,log=TRUE)-log(2)
	return(pdf)
}

pdweibull=function (x, c=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE&x<=mu]=0.5*exp(-(mu-x[x<=mu])**c/sigma**c)
        cdf[log.p==FALSE&lower.tail==TRUE&x>mu]=1-0.5*exp(-(x[x>mu]-mu)**c/sigma**c)
        cdf[log.p==FALSE&lower.tail==FALSE&x<=mu]=1-0.5*exp(-(mu-x[x<=mu])**c/sigma**c)
        cdf[log.p==FALSE&lower.tail==FALSE&x>mu]=0.5*exp(-(x[x>mu]-mu)**c/sigma**c)
        cdf[log.p==TRUE&lower.tail==TRUE&x<=mu]=-(mu-x[x<=mu])**c/sigma**c-log(2)
        cdf[log.p==TRUE&lower.tail==TRUE&x>mu]=log(1-0.5*exp(-(x[x>mu]-mu)**c/sigma**c))
        cdf[log.p==TRUE&lower.tail==FALSE&x<=mu]=log(1-0.5*exp(-(mu-x[x<=mu])**c/sigma**c))
        cdf[log.p==TRUE&lower.tail==FALSE&x>mu]=-(x[x>mu]-mu)**c/sigma**c-log(2)
	return(cdf)
}

vardweibull=function (p, c=1, mu=0, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
        var=p
	var[p<=0.5]=mu-sigma*(-log(2*p[p<=0.5]))**(1/c)
	var[p>0.5]=mu+sigma*(-log(2*(1-p[p>0.5])))**(1/c)
	return(var)
}

esdweibull=function (p, c=1, mu=0, sigma=1)
{
	f=function (x) {vardweibull(x,c=c,mu=mu,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Exponentiated Weibull

dexpweibull=function (x, a=1, alpha=1, sigma=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=a*dweibull(x,shape=alpha,scale=sigma)*(1-pweibull(x,shape=alpha,scale=sigma))**(a-1)
        pdf[log==TRUE]=log(a)+dweibull(x,shape=alpha,scale=sigma,log=TRUE)+(a-1)*log(1-pweibull(x,shape=alpha,scale=sigma))
	return(pdf)
}

pexpweibull=function (x, a=1, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(1-pweibull(x,shape=alpha,scale=sigma))**a
        cdf[log.p==FALSE&lower.tail==FALSE]=1-(1-pweibull(x,shape=alpha,scale=sigma))**a
	cdf[log.p==TRUE&lower.tail==TRUE]=a*log(1-pweibull(x,shape=alpha,scale=sigma))
	cdf[log.p==TRUE&lower.tail==FALSE]=log(1-(1-pweibull(x,shape=alpha,scale=sigma))**a)
	return(cdf)
}

varexpweibull=function (p, a=1, alpha=1, sigma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=sigma*(-log(1-p**(1/a)))**(1/alpha)
	return(var)
}

esexpweibull=function (p, a=1, alpha=1, sigma=1)
{
	f=function (x) {varexpweibull(x,a=a,alpha=alpha,sigma=sigma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Generalized power Weibull

dgenpowerweibull=function (x, a=1, theta=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=a*theta*x**(a-1)*(1+x**a)**(theta-1)*exp(1-(1+x**a)**theta)
        pdf[log==TRUE]=log(a*theta)+(a-1)*log(x)+(theta-1)*log(1+x**a)+1-(1+x**a)**theta
	return(pdf)
}

pgenpowerweibull=function (x, a=1, theta=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-exp(1-(1+x**a)**theta)
        cdf[log.p==FALSE&lower.tail==FALSE]=exp(1-(1+x**a)**theta)
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1-exp(1-(1+x**a)**theta))
	cdf[log.p==TRUE&lower.tail==FALSE]=1-(1+x**a)**theta
	return(cdf)
}

vargenpowerweibull=function (p, a=1, theta=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=((1-log(1-p))**(1/theta)-1)**(1/a)
	return(var)
}

esgenpowerweibull=function (p, a=1, theta=1)
{
	f=function (x) {vargenpowerweibull(x,a=a,theta=theta)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Chen

dchen=function (x, b=1, lambda=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=b*lambda*x**(b-1)*exp(x**b)*exp(lambda-lambda*exp(x**b))
        pdf[log==TRUE]=log(b*lambda)+(b-1)*log(x)+x**b+lambda-lambda*exp(x**b)
	return(pdf)
}

pchen=function (x, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-exp(lambda-lambda*exp(x**b))
        cdf[log.p==FALSE&lower.tail==FALSE]=exp(lambda-lambda*exp(x**b))
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1-exp(lambda-lambda*exp(x**b)))
	cdf[log.p==TRUE&lower.tail==FALSE]=lambda-lambda*exp(x**b)
	return(cdf)
}

varchen=function (p, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(log(1-(1/lambda)*log(1-p)))**(1/b)
	return(var)
}

eschen=function (p, b=1, lambda=1)
{
	f=function (x) {varchen(x,b=b,lambda=lambda)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Xie


dxie=function (x, a=1, b=1, lambda=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=b*lambda*(x/a)**(b-1)*exp(a*lambda+(x/a)**b-a*lambda*exp((x/a)**b))
        pdf[log==TRUE]=log(b*lambda)+(b-1)*log(x/a)+a*lambda+(x/a)**b-a*lambda*exp((x/a)**b)
	return(pdf)
}

pxie=function (x, a=1, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-exp(a*lambda-a*lambda*exp((x/a)**b))
        cdf[log.p==FALSE&lower.tail==FALSE]=exp(a*lambda-a*lambda*exp((x/a)**b))
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1-exp(a*lambda-a*lambda*exp((x/a)**b)))
	cdf[log.p==TRUE&lower.tail==FALSE]=a*lambda-a*lambda*exp((x/a)**b)
	return(cdf)
}

varxie=function (p, a=1, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=a*(log(1-(log(1-p))/(lambda*a)))**(1/b)
	return(var)
}

esxie=function (p, a=1, b=1, lambda=1)
{
	f=function (x) {varxie(x,a=a,b=b,lambda=lambda)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Tukey Lambda

varTL=function (p, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(p**lambda-(1-p)**lambda)/lambda
	return(var)
}

esTL=function (p, lambda=1)
{
	f=function (x) {varTL(x,lambda=lambda)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Ramberg Schmeiser

varRS=function (p, b=1, c=1, d=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(p**b-(1-p)**c)/d
	return(var)
}

esRS=function (p, b=1, c=1, d=1)
{
	f=function (x) {varRS(x,b=b,c=c,d=d)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Freimer


varFR=function (p, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=(1/a)*((p**b-1)/b-((1-p)**c-1)/c)
	return(var)
}

esFR=function (p, a=1, b=1, c=1)
{
	f=function (x) {varFR(x,a=a,b=b,c=c)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Hankin Lee

varHL=function (p, a=1, b=1, c=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=c*p**a/(1-p)**b
	return(var)
}

esHL=function (p, a=1, b=1, c=1)
{
	f=function (x) {varHL(x,a=a,b=b,c=c)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Loglog

dloglog=function (x, a=1, lambda=2, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=a*log(lambda)*x**(a-1)*lambda**(x**a)*exp(1-lambda**(x**a))
        pdf[log==TRUE]=log(a)+log(log(lambda))+(a-1)*log(x)+(x**a)*log(lambda)+1-lambda**(x**a)
	return(pdf)
}

ploglog=function (x, a=1, lambda=2, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-exp(1-lambda**(x**a))
        cdf[log.p==FALSE&lower.tail==FALSE]=exp(1-lambda**(x**a))
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1-exp(1-lambda**(x**a)))
	cdf[log.p==TRUE&lower.tail==FALSE]=1-lambda**(x**a)
	return(cdf)
}


varloglog=function (p, a=1, lambda=2, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=((log(1-log(1-p)))/log(lambda))**(1/a)
	return(var)
}

esloglog=function (p, a=1, lambda=2)
{
	f=function (x) {varloglog(x,a=a,lambda=lambda)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Exponential logarithmic

dexplog=function (x, a=0.5, b=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=-b*(1-a)*(exp(-b*x))/(log(a)*(1-(1-a)*exp(-b*x)))
        pdf[log==TRUE]=log(b)+log(1-a)-b*x-log(-log(a))-log(1-(1-a)*exp(-b*x))
	return(pdf)
}

pexplog=function (x, a=0.5, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=1-(log(1-(1-a)*exp(-b*x)))/log(a)
        cdf[log.p==FALSE&lower.tail==FALSE]=(log(1-(1-a)*exp(-b*x)))/log(a)
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1-(log(1-(1-a)*exp(-b*x)))/log(a))
	cdf[log.p==TRUE&lower.tail==FALSE]=log(-log(1-(1-a)*exp(-b*x)))-log(-log(a))
	return(cdf)
}

varexplog=function (p, a=0.5, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=-(1/b)*log((1-a**(1-p))/(1-a))
	return(var)
}

esexplog=function (p, a=0.5, b=1)
{
	f=function (x) {varexplog(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Exponential geometric

dexpgeo=function (x, theta=0.5, lambda=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=lambda*theta*exp(-lambda*x)*(1-(1-theta)*exp(-lambda*x))**(-2)
        pdf[log==TRUE]=log(lambda*theta)-lambda*x-2*log(1-(1-theta)*exp(-lambda*x))
	return(pdf)
}

pexpgeo=function (x, theta=0.5, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=theta*exp(-lambda*x)*(1-(1-theta)*exp(-lambda*x))**(-1)
        cdf[log.p==FALSE&lower.tail==FALSE]=(1-exp(-lambda*x))*(1-(1-theta)*exp(-lambda*x))**(-1)
	cdf[log.p==TRUE&lower.tail==TRUE]=log(theta)-lambda*x-log(1-(1-theta)*exp(-lambda*x))
	cdf[log.p==TRUE&lower.tail==FALSE]=log(1-exp(-lambda*x))-log(1-(1-theta)*exp(-lambda*x))
	return(cdf)
}

varexpgeo=function (p, theta=0.5, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=-(1/lambda)*log(p/(theta+(1-theta)*p))
	return(var)
}

esexpgeo=function (p, theta=0.5, lambda=1)
{
	f=function (x) {varexpgeo(x,theta=theta,lambda=lambda)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Exponential Poisson

dexppois=function (x, b=1, lambda=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=b*lambda*(exp(-b*x-lambda+lambda*exp(-b*x)))/(1-exp(-lambda))
        pdf[log==TRUE]=log(b*lambda)-b*x-lambda+lambda*exp(-b*x)-log(1-exp(-lambda))
	return(pdf)
}

pexppois=function (x, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(1-exp(-lambda+lambda*exp(-b*x)))/(1-exp(-lambda))
        cdf[log.p==FALSE&lower.tail==FALSE]=(exp(-lambda+lambda*exp(-b*x))-exp(-lambda))/(1-exp(-lambda))
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1-exp(-lambda+lambda*exp(-b*x)))-log(1-exp(-lambda))
	cdf[log.p==TRUE&lower.tail==FALSE]=log(exp(-lambda+lambda*exp(-b*x))-exp(-lambda))-log(1-exp(-lambda))
	return(cdf)
}

varexppois=function (p, b=1, lambda=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=-(1/b)*log((1/lambda)*log(1-p+p*exp(-lambda))+1)
	return(var)
}

esexppois=function (p, b=1, lambda=1)
{
	f=function (x) {varexppois(x,b=b,lambda=lambda)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Topp Leone

dTL2=function (x, b=1, log=FALSE)
{
	pdf=x
        pdf[log==FALSE]=2*b*(x*(2-x))**(b-1)*(1-x)
        pdf[log==TRUE]=log(2*b)+(b-1)*log(x*(2-x))+log(1-x)
	return(pdf)
}

pTL2=function (x, b=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(x*(2-x))**b
        cdf[log.p==FALSE&lower.tail==FALSE]=1-(x*(2-x))**b
	cdf[log.p==TRUE&lower.tail==TRUE]=b*log(x*(2-x))
	cdf[log.p==TRUE&lower.tail==FALSE]=log(1-(x*(2-x))**b)
	return(cdf)
}

varTL2=function (p, b=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=1-sqrt(1-p**(1/b))
	return(var)
}

esTL2=function (p, b=1)
{
	f=function (x) {varTL2(x,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Quadratic

dquad=function (x, a=0, b=1, log=FALSE)
{
	alpha=12/(b-a)**3
        beta=(a+b)/2
        pdf=x
        pdf[log==FALSE]=alpha*(x-beta)**2
        pdf[log==TRUE]=log(alpha)+2*log(x-beta)
	return(pdf)
}

pquad=function (x, a=0, b=1, log.p=FALSE, lower.tail=TRUE)
{
        alpha=12/(b-a)**3
        beta=(a+b)/2
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(alpha/3)*((x-beta)**3+(beta-a)**3)
        cdf[log.p==FALSE&lower.tail==FALSE]=1-(alpha/3)*((x-beta)**3+(beta-a)**3)
	cdf[log.p==TRUE&lower.tail==TRUE]=log(alpha/3)+log((x-beta)**3+(beta-a)**3)
	cdf[log.p==TRUE&lower.tail==FALSE]=log(1-(alpha/3)*((x-beta)**3+(beta-a)**3))
	return(cdf)
}

varquad=function (p, a=0, b=1, log.p=FALSE, lower.tail=TRUE)
{
        alpha=12/(b-a)**3
        beta=(a+b)/2
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=beta+sign(3*p/alpha-(beta-a)**3)*(abs(3*p/alpha-(beta-a)**3))**(1/3)
	return(var)
}

esquad=function (p, a=0, b=1)
{
	f=function (x) {varquad(x,a=a,b=b)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Schabe

dschabe=function (x, gamma=0.5, theta=1, log=FALSE)
{
        pdf=x
        pdf[log==FALSE]=(1/theta)*(2*gamma+(1-gamma)*x/theta)*(gamma+x/theta)**(-2)
        pdf[log==TRUE]=-log(theta)+log(2*gamma+(1-gamma)*x/theta)-2*log(gamma+x/theta)
	return(pdf)
}

pschabe=function (x, gamma=0.5, theta=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=(1+gamma)*x/(x+gamma*theta)
        cdf[log.p==FALSE&lower.tail==FALSE]=gamma*(theta-x)/(x+gamma*theta)
	cdf[log.p==TRUE&lower.tail==TRUE]=log(1+gamma)+log(x)-log(x+gamma*theta)
	cdf[log.p==TRUE&lower.tail==FALSE]=log(gamma)+log(theta-x)-log(x+gamma*theta)
	return(cdf)
}

varschabe=function (p, gamma=0.5, theta=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=p*gamma*theta/(1+gamma-p)
	return(var)
}

esschabe=function (p, gamma=0.5, theta=1)
{
	f=function (x) {varschabe(x,gamma=gamma,theta=theta)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Birnbaum Saunders

dBS=function (x, gamma=1, log=FALSE)
{
        pdf=x
        pdf[log==FALSE]=(x**(1/2)+x**(-1/2))*(dnorm((x**(1/2)-x**(-1/2))/gamma))/(2*gamma*x)
        pdf[log==TRUE]=log(x**(1/2)+x**(-1/2))+dnorm((x**(1/2)-x**(-1/2))/gamma,log=TRUE)-log(2*gamma)-log(x)
	return(pdf)
}

pBS=function (x, gamma=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=pnorm((x**(1/2)-x**(-1/2))/gamma,log.p=log.p,lower.tail=lower.tail)
	return(cdf)
}

varBS=function (p, gamma=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=0.25*(gamma*qnorm(p)+sqrt(4+gamma**2*(qnorm(p))**2))**2
	return(var)
}

esBS=function (p, gamma=1)
{
	f=function (x) {varBS(x,gamma=gamma)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}

#Generalized extreme value

dgev=function (x, mu=0, sigma=1, xi=1, log=FALSE)
{
        pdf=x
        pdf[log==FALSE]=(1/sigma)*(1+xi*(x-mu)/sigma)**(-1/xi-1)*exp(-(1+xi*(x-mu)/sigma)**(-1/xi))
        pdf[log==TRUE]=-log(sigma)-(1/xi+1)*log(1+xi*(x-mu)/sigma)-(1+xi*(x-mu)/sigma)**(-1/xi)
	return(pdf)
}

pgev=function (x, mu=0, sigma=1, xi=1, log.p=FALSE, lower.tail=TRUE)
{
	cdf=x
        cdf[log.p==FALSE&lower.tail==TRUE]=exp(-(1+xi*(x-mu)/sigma)**(-1/xi))
        cdf[log.p==FALSE&lower.tail==FALSE]=1-exp(-(1+xi*(x-mu)/sigma)**(-1/xi))
	cdf[log.p==TRUE&lower.tail==TRUE]=-(1+xi*(x-mu)/sigma)**(-1/xi)
	cdf[log.p==TRUE&lower.tail==FALSE]=log(1-exp(-(1+xi*(x-mu)/sigma)**(-1/xi)))
	return(cdf)
}

vargev=function (p, mu=0, sigma=1, xi=1, log.p=FALSE, lower.tail=TRUE)
{
	if (log.p==TRUE) p=exp(p)
        if (lower.tail==FALSE) p=1-p
	var=mu-(sigma/xi)+(sigma/xi)*(-log(p))**(-xi)
	return(var)
}

esgev=function (p, mu=0, sigma=1, xi=1)
{
	f=function (x) {vargev(x,mu=mu,sigma=sigma,xi=xi)}
        es=p
        for (i in 1:length(p)) {es[i]=(1/p[i])*integrate(f,lower=0,upper=p[i],stop.on.error=FALSE)$value}
	return(es)
}
