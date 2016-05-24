

#Gompertz probability density function

dgompertz=function(x,alpha=1,beta=1)
{
	ret = ifelse(x<=0 | alpha<=0 | beta<=0 ,
          NaN,alpha*exp(beta*x)*exp((1/beta)*alpha*(1-exp(beta*x))))
	return(ret)
}


#Makeham probability density function

dmakeham=function(x,alpha=1,beta=1,epsilon=1)
{
	ret = ifelse(x<=0 | alpha<=0 | beta<=0 | epsilon<=0,NaN,
          epsilon+alpha*exp(beta*x)*
          exp(-x*epsilon+(1/beta)*alpha*(1-exp(beta*x))))
	return(ret)
}



#Perks probability density function

dperks=function(x,alpha=1,beta=1)
{
	ret	= ifelse(x<=0 | alpha<=0 | beta<=0,NaN,
           alpha*exp(beta*x)*(1+alpha)/
           (1+alpha*exp(beta*x))**2)
	return(ret)
}



#Beard probability density function

dbeard=function(x,alpha=1,beta=1,rho=1)
{
	ret=ifelse(x<=0 | alpha<=0 | beta<=0 | rho<=0,NaN,
           alpha*exp(beta*x)*(1+alpha*rho)**(rho**(-1/beta))/
           (1+alpha*rho*exp(beta*x))**(1+rho**(-1/beta)))
	return(ret)
}



#Makeham-Perks probability density function

dmakehamperks=function(x,alpha=1,beta=1,epsilon=1)
{
	ret = ifelse(x<=0 | alpha<=0 | beta<=0 | epsilon<=0,NaN,
          (epsilon+alpha*exp(beta*x))/(1+alpha*exp(beta*x))*
          exp(-epsilon*x-(1/beta)*(epsilon-1)*log((1+alpha)/(1+alpha*exp(beta*x)))))
	return(ret)
}



#Makeham-Beard probability density function

dmakehambeard=function(x,alpha=1,beta=1,rho=1,epsilon=1)
{
	ret = ifelse(x<=0 | alpha<=0 | beta<=0 | epsilon<=0 | rho<=0,NaN,
          (epsilon+alpha*exp(beta*x))/(1+alpha*rho*exp(beta*x))*
          exp(-epsilon*x-(1/beta)*(1/rho)*(rho*epsilon-1)*log((1+alpha*rho)/(1+alpha*rho*exp(beta*x)))))
	return(ret)
}



#Exponential probability density function
#
#dexponential=function(x,alpha=1)
#{
#	ret	= ifelse(x<=0 | alpha<=0, NaN,alpha*exp(-alpha*x))
#	return(ret)
#}



#Pareto probability density function

dpareto=function(x,alpha=1,m=1)
{
	ret	= ifelse(x<=m | m<=0 | alpha<=0,NaN,(m/x)**(alpha)*alpha/x)
	return(ret)
}



#Weibull probability density function
#
#dweibull=function(x,alpha=1,sigma=1)
#{
#	ret	= ifelse(x<=0 | sigma<=0 | alpha<=0,NaN,
#           sigma*alpha**(-sigma)*x**(sigma-1)*exp(-(x/alpha)**sigma))
#	return(ret)
#}



#Logistic probability density function

#dlogistic=function(x,alpha=1,sigma=1)
#{
#	ret	= ifelse (sigma<=0,NaN,dlogis(x,location=alpha,scale=sigma))
#	return(ret)
#}
#
#The function {\sf dlogis} is  from the {\sf R} base package.



#Log-logistic probability density function

dloglogistic=function(x,alpha=1,sigma=1)
{
	ret	= ifelse(x<=0 | alpha<=0 | sigma<=0, NaN,alpha*sigma*x**(sigma-1)*(1+alpha*x**sigma)**(-2))
	return(ret)
}



#Normal probability density function
#
#dnormal=function(x,alpha=1,sigma=1)
#{
#	ret	= ifelse (sigma<=0,NaN,dnorm(x,mean=alpha,sd=sigma))
#	return(ret)
#}
#
#The function {\sf dnorm} is  from the {\sf R} base package.


#Lognormal probability density function
#
#dlognormal=function(x,alpha=1,sigma=1)
#{
#	ret	= ifelse(x<=0 | sigma<=0, NaN,
#           dlnorm(x,meanlog=alpha,sdlog=sigma))
#	return(ret)
#}
#
#The function {\sf dlnorm} is  from the {\sf R} base package.


#Inverse-Gaussian probability density function

dinvgauss=function(x,alpha=1,sigma=1)
{
	ret	= ifelse(x<=0 | alpha<=0 | sigma<=0, NaN,
           sqrt(sigma/(2*pi*x**3))*exp(-sigma*(x-alpha)**2/(2*alpha**2*x)))
	return(ret)
}



#Gamma probability density function
#
#dgammad=function(x,alpha=1,lambda=1)
#{
#	ret	= ifelse(x<=0 | alpha<=0 | lambda<=0, NaN,
#           dgamma(x,shape=lambda,scale=alpha))
#	return(ret)
#}
#
#The function {\sf dgamma} is  from the {\sf R} base package.



#Generalized-gamma probability density function

dgengammad=function(x,b=1,d=1,k=1)
{
	ret	= ifelse(x<=0 | b<=0 | d<=0 | k<=0, NaN,
            d*(b**(-d*k))*x**(d*k-1)*exp(-(x/b)**d)/gamma(k))           
	return(ret)
}


#Linear probability density function

dlinear=function(x,a=1,b=1)
{
	ret	= ifelse(x<=0 | a<0 | b<0 | a==0&b==0, NaN, (a+b*x)*exp(-a*x-b*x*x/2))
	return(ret)
}




#Uniform probability density function
#
#duniform=function(x,a=1,b=2)
#{
#	ret	= ifelse(b<=a | a<=0,NaN,1/(b-a))
#	return(ret)
#}




#Kumaraswamy probability density function

dkum=function(x,a=1,b=1)
{
	ret	= ifelse(x<=0 | x>=1 | a<=0 | b<=0,
           NaN,a*b*x**(a-1)*(1-x**a)**(b-1))
	return(ret)
}



#Gumbel II probability density function

dfrechet=function(x,a=1,b=1)
{
	ret	= ifelse(x<=0 | a<=0 | b<=0,
           NaN, a*b*x**(a-1)*exp(-b*x**(-a)))
	return(ret)
}



#Hjorth probability density function

dhjorth=function(x,delta=1,theta=1,beta=1)
{
	ret	= ifelse(x<=0 | delta<=0 | theta<=0 | beta<=0,
           NaN, ((1+beta*x)*delta*x+theta)*(1+beta*x)**(-1-theta/beta)*exp(-delta*x*x/2))
	return(ret)
}



#Additive Weibull probability density function

daddweibull=function(x,a=1,b=1,c=1,d=1)
{
	ret	= ifelse(x<=0 | a<=0 | b<=0 | c<=0 | d<=0,
           NaN, (a**b*b*x**{b-1}+c**d*d*x**{d-1})*exp(-(a*x)**b-(c*x)**d))
	return(ret)
}



#Schabe probability density function

dschabe=function(x,theta=1,gamma=0.5)
{
	ret	= ifelse(x>theta | gamma<=0 | gamma>=1 | theta<=0,
           NaN,(2*gamma+(1-gamma)*x/theta)/(theta*(gamma+x/theta)**2))
	return(ret)
}




#Lai probability density function

dlai=function(x,lambda=1,beta=1,nu=1)
{
	ret	= ifelse(x<=0 | lambda<=0 | beta<0 | nu<0,
           NaN,lambda*(beta+nu*x)*x**(beta-1)*exp(nu*x)*
           exp(-lambda*x**beta*exp(nu*x)))
	return(ret)
}




#Xie probability density function

dxie=function(x,lambda=1,alpha=1,beta=1)
{
	ret	= ifelse(x<=0 | lambda<=0 | alpha<=0 | beta<=0,
           NaN,
           lambda*beta*(x/alpha)**(beta-1)*
           exp((x/alpha)**beta+lambda*alpha*(1-exp((x/alpha)**beta))))
	return(ret)
}




#$J$-shaped probability density function

djshape=function(x,b=1,nu=1)
{
	y=1-x/b
    ret	= ifelse(y<=0 | y>=1 | b<=0 | nu<=0,
           NaN,(2/b)*nu*y*(1-y*y)**(nu-1))
	return(ret)
}




#Beta probability density function
#
#dbetad=function(x,a=1,b=1)
#{
#    ret	= ifelse(x<=0 | x>=1 | a<=0 | b<=0,
#           NaN,dbeta(x,shape1=a,shape2=b))
#	return(ret)
#}
#
#The function  {\sf dbeta} is  from the {\sf R} base package.



#Exponentiated exponential probability density function

dee=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0,
           NaN,dgen.exp(x,
           alpha=alpha,lambda=lambda))
	return(ret)
}

#The function {\sf dgen.exp} is  from the {\sf R} contributed package {\sf reliaR}.



#Inverse generalized exponential probability density function

dige=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0,
           NaN,dinv.genexp(x,
           alpha=alpha,lambda=lambda))
	return(ret)
}

#The function {\sf dinv.genexp} is  from the {\sf R} contributed package {\sf reliaR}.



#Exponentiated Weibull probability density function

dew=function(x,alpha=1,c=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | c<=0 | lambda<=0,
           NaN,lambda*dexpo.weibull(lambda*x,
           alpha=c,theta=alpha))
	return(ret)
}

#The function {\sf dexpo.weibull} is  from the {\sf R} contributed package {\sf reliaR}.



#Exponentiated logistic probability density function

del=function(x,alpha=1,beta=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | beta<=0, NaN,
           dexpo.logistic(x,
           alpha=alpha,beta=beta))
	return(ret)
}

#The function {\sf dexpo.logistic} is  from the {\sf R} contributed package {\sf reliaR}.



#BurrX probability density function

dburrx=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0, NaN,
           dburrX(x,alpha=alpha,
           lambda=lambda))
	return(ret)
}

#The function {\sf dburrX} is  from the {\sf R} contributed package {\sf reliaR}.



#Gumbel probability density function

dgumbeld=function(x,mu=1,sigma=1)
{
    ret	= ifelse(sigma<=0, NaN,
           dgumbel(x,mu=mu,sigma=sigma))
	return(ret)
}

#The function {\sf dgumbel} is  from the {\sf R} contributed package {\sf reliaR}.



#Exponential extension probability density function

dexpext=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0,
           NaN, alpha*lambda*(1+lambda*x)**(alpha-1)*
           exp(1-(1+lambda*x)**alpha))
	return(ret)
}




#Chen probability density function

dchen=function(x,beta=1,lambda=1)
{
    ret	= ifelse(x<=0 | beta<=0 | lambda<=0,
           NaN, beta*lambda*
           x**(beta-1)*exp(x**beta)*exp(lambda-lambda*exp(x**beta)))
	return(ret)
}





#Loggamma probability density function

dlgammad=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=1 | alpha<=0 | lambda<=0, NaN,
           dlgamma(x,shapelog=alpha,
           ratelog=lambda))
	return(ret)
}

#The function {\sf dlgamma} is  from the {\sf R} contributed package {\sf actuar}.



#Marshall-Olkin exponential probability density function

dmoe=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0,
           NaN, lambda*
           exp(lambda*x)/
           (exp(lambda*x)-1+alpha)**2)
	return(ret)
}




#Marshall-Olkin Weibull probability density function

dmow=function(x,alpha=1,beta=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | beta<=0 | lambda<=0,
           NaN, beta*lambda**beta*x**{beta-1}*
           exp((lambda*x)**beta)/
           (exp((lambda*x)**beta)-1+alpha)**2)
	return(ret)
}




#Logistic exponential probability density function

dle=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0, NaN,
           dlogis.exp(x,
           alpha=alpha,lambda=lambda))
	return(ret)
}

#The function {\sf dlogis.exp} is  from the {\sf R} contributed package {\sf reliaR}.



#Logistic Rayleigh probability density function

dlr=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0, NaN,
           dlogis.rayleigh(x,
           alpha=alpha,lambda=lambda))
	return(ret)
}

#The function {\sf dlogis.rayleigh} is  from the {\sf R} contributed package {\sf reliaR}.



#Loglog probability density function

dlld=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0,
           NaN,alpha*log(lambda)*
           x**(alpha-1)*lambda**(x**alpha)*exp(1-lambda**(x**alpha)))
	return(ret)
}




#Generalized power Weibull probability density function

dgpw=function(x,alpha=1,theta=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | theta<=0,
           NaN,alpha*theta*x**{alpha-1}*
           (1+x**alpha)**{theta-1}*exp(1-(1+x**alpha)**theta))
	return(ret)
}




#Flexible Weibull probability density function

dfw=function(x,alpha=1,beta=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | beta<=0,
           NaN,(alpha+beta/(x*x))*
           exp(alpha*x-beta/x)*exp(-exp(alpha*x-beta/x)))
	return(ret)
}




#Generalized $F$ probability density function

dgenF=function(x,beta=0,sigma=1,m1=1,m2=1)
{
    ret	= ifelse(x<=0 | sigma<=0 | m1<=0 | m2<=0,
           NaN,exp(-beta*m1/sigma)*(m1/m2)**m1*x**(m1/sigma-1)/(sigma*beta(m1,m2)*(1+(m1/m2)*(exp(-beta)*x)**(1/sigma))**(m1+m2)))
	return(ret)
}



#Gompertz hazard rate function

hgompertz=function(x,alpha=1,beta=1)
{
	ret	= ifelse(x<=0 | alpha<=0 | beta<=0,
           NaN,alpha*exp(beta*x))
	return(ret)
}



#Makeham hazard rate function

hmakeham=function(x,alpha=1,beta=1,epsilon=1)
{
	ret	=ifelse(x<=0 | alpha<=0 | beta<=0 | epsilon <=0,NaN,
          epsilon+alpha*exp(beta*x))
	return(ret)
}



#Perks hazard rate function

hperks=function(x,alpha=1,beta=1)
{
	ret	= ifelse(x<=0 | alpha<=0 | beta<=0,NaN,
           alpha*exp(beta*x) /
           (1+alpha*exp(beta*x)))
	return(ret)
}



#Beard hazard rate function

hbeard=function(x,alpha=1,beta=1,rho=1)
{
	ret=ifelse(x<=0 | alpha<=0 | beta<=0 | rho<=0,NaN,
           alpha*exp(beta*x) /
           (1+alpha*rho*exp(beta*x)))
	return(ret)
}



#Makeham-Perks hazard rate function

hmakehamperks=function(x,alpha=1,beta=1,epsilon=1)
{
	ret	= ifelse(x<=0 | alpha<=0 | beta<=0 | epsilon<=0,
           NaN,(epsilon+alpha*exp(beta*x)) /
           (1+alpha*exp(beta*x)))
	return(ret)
}



#Makeham-Beard hazard rate function

hmakehambeard=function(x,alpha=1,beta=1,rho=1,epsilon=1)
{
	ret	= ifelse(x<=0 | alpha<=0 | beta<=0 | rho<=0 | epsilon<=0,
           NaN,(epsilon+alpha*exp(beta*x)) /
           (1+alpha*rho*exp(beta*x)))
	return(ret)
}



#Exponential hazard rate function

hexponential=function(x,alpha=1)
{
	ret	= ifelse(x<=0 | alpha<=0, NaN,alpha)
	return(ret)
}



#Pareto hazard rate function

hpareto=function(x,alpha=1,m=1)
{
	ret	= ifelse(x<=m | m<=0 | alpha<=0,NaN,alpha/x)
	return(ret)
}



#Weibull hazard rate function

hweibull=function(x,alpha=1,sigma=1)
{
	ret	= ifelse(x<=0 | sigma<=0 | alpha<=0,NaN,
           alpha*x**(alpha-1)/sigma**alpha)
	return(ret)
}



#Logistic hazard rate function

hlogistic=function(x,alpha=1,sigma=1)
{
	ret	= ifelse(sigma<=0,NaN,1 / (sigma*
           (1+exp( -(x-alpha)/sigma))))
	return(ret)
}



#Log-logistic hazard rate function

hloglogistic=function(x,alpha=1,sigma=1)
{
	ret	= ifelse(x<=0 | sigma<=0 | alpha<=0, NaN,
                alpha*sigma*x**(sigma-1)/(1+alpha*x**sigma))
	return(ret)
}



#Normal hazard rate function

hnormal=function(x,alpha=1,sigma=1)
{
	ret	= ifelse(x<=0 | sigma<=0,NaN,dnorm(x,mean=alpha,sd=sigma) /
          (1-pnorm(x,mean=alpha,sd=sigma)))
	return(ret)
}

#The functions {\sf dnorm} and {\sf pnorm} are  from the {\sf R} base package.


#Lognormal hazard rate function

hlognormal=function(x,alpha=1,sigma=1)
{
	ret	= ifelse(x<=0 | sigma<=0, NaN,
           dlnorm(x,meanlog=alpha,sdlog=sigma) /
           (1-plnorm(x,meanlog=alpha,sdlog=sigma)))
	return(ret)
}

#The functions {\sf dlnorm} and {\sf plnorm} are  from the {\sf R} base package.


#Inverse-Gaussian hazard rate function

hinvgauss=function(x,alpha=1,sigma=1)
{
	ret=NaN
        if (x>0&alpha>0&sigma>0)
        {tt=sqrt(sigma/(2*pi*x**3))*exp(-sigma*(x-alpha)**2/(2*alpha**2*x))
        tt2=pnorm(sqrt(sigma/x)*(x/alpha-1))+exp(2*sigma/alpha)*pnorm(-sqrt(sigma/x)*(x/alpha+1))
        ret=tt/(1-tt2)}
	return(ret)
}

#The function {\sf pnorm} is from the {\sf R} base package.


#Gamma hazard rate function

hgamma=function(x,alpha=1,lambda=1)
{
	ret	= ifelse(x<=0 | alpha<=0 | lambda<=0, NaN,
           dgamma(x,shape=lambda,scale=alpha) /
           (1-pgamma(x,shape=lambda,scale=alpha)))
	return(ret)
}

#The functions {\sf dgamma} and {\sf pgamma} are  from the {\sf R} base package.



#Generalized-gamma hazard rate function

hgengamma=function(x,b=1,d=1,k=1)
{
	ret	= ifelse(x<=0 | b<=0 | d<=0 | k<=0, NaN,
          d*(b**(-d*k))*x**(d*k-1)*exp(-(x/b)**d)/(gamma(k)*(1-pgamma((x/b)**d,shape=k))))
	return(ret)
}

#The function {\sf pgamma} is from the {\sf R} base package.


#Linear hazard rate function

hlinear=function(x,a=1,b=1)
{
	ret	= ifelse(x<=0 | a<0 | b<0 | a==0&b==0, NaN, a+b*x)
	return(ret)
}




#Uniform hazard rate function

huniform=function(x,a=1,b=2)
{
	ret	= ifelse(x<=a | x>=b | b<=a | a<=0,
           NaN, 1/(b-x))
	return(ret)
}




#Kumaraswamy hazard rate function

hkum=function(x,a=1,b=1)
{
	ret	= ifelse(x<=0 | x>=1 | a<=0 | b<=0,
           NaN,a*b*x**(a-1)/(1-x**a))
	return(ret)
}



#Gumbel II hazard rate function

hfrechet=function(x,a=1,b=1)
{
	ret	= ifelse(x<=0 | a<=0 | b<=0,
           NaN, a*b*x**(a-1))
	return(ret)
}



#Hjorth hazard rate function

hhjorth=function(x,delta=1,theta=1,beta=1)
{
	ret	= ifelse(x<=0 | delta<=0 | theta<=0 | beta<=0,
           NaN, delta*x+theta/(1+beta*x))
	return(ret)
}



#Additive Weibull hazard rate function

haddweibull=function(x,a=1,b=1,c=1,d=1)
{
	ret	= ifelse(x<=0 | a<=0 | b<=0 | c<=0 | d<=0,
           NaN, a**b*b*x**{b-1}+
           c**d*d*x**{d-1})
	return(ret)
}



#Schabe hazard rate function

hschabe=function(x,theta=1,gamma=0.5)
{
	ret	= ifelse(x>theta | gamma<=0 | gamma>=1 | theta<=0,
           NaN,1/(theta*gamma+x)+
           1/(theta-x))
	return(ret)
}




#Lai hazard rate function

hlai=function(x,lambda=1,beta=1,nu=1)
{
	ret	= ifelse(x<=0 | lambda<=0 | beta<0 | nu<0,
           NaN,lambda*(beta+nu*x)*
           x**(beta-1)*exp(nu*x))
	return(ret)
}




#Xie hazard rate function

hxie=function(x,lambda=1,alpha=1,beta=1)
{
	ret	= ifelse(x<=0 | lambda<=0 | alpha<=0 | beta<=0,
           NaN,
           lambda*beta*(x/alpha)**(beta-1)*exp((x/alpha)**beta))
	return(ret)
}




#$J$-shaped hazard rate function

hjshape=function(x,b=1,nu=1)
{
	y=1-x/b
    ret	= ifelse(y<=0 | y>=1 | b<=0 | nu<=0,
           NaN,(2/b)*nu*y*(1-y*y)**(nu-1)/
           (1-(1-y*y)**nu))
	return(ret)
}




#Beta hazard rate function

hbeta=function(x,a=1,b=1)
{
    ret	= ifelse(x<=0 | x>=1 | a<=0 | b<=0,
           NaN,dbeta(x,shape1=a,shape2=b)/
           (1-pbeta(x,shape1=a,shape2=b)))
	return(ret)
}

#The functions  {\sf dbeta} and {\sf pbeta} are  from the {\sf R} base package.



#Exponentiated exponential hazard rate function

hee=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0,
           NaN,dgen.exp(x,
           alpha=alpha,lambda=lambda)/
           (1-pgen.exp(x,
           alpha=alpha,lambda=lambda)))
	return(ret)
}

#The functions {\sf dgen.exp} and {\sf pgen.exp} are  from the {\sf R} contributed package {\sf reliaR}.



#Inverse generalized exponential hazard rate function

hige=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0,
           NaN,dinv.genexp(x,
           alpha=alpha,lambda=lambda)/
           (1-pinv.genexp(x,
           alpha=alpha,lambda=lambda)))
	return(ret)
}

#The functions {\sf dinv.genexp} and {\sf pinv.genexp} are  from the {\sf R} contributed package {\sf reliaR}.



#Exponentiated Weibull hazard rate function

hew=function(x,alpha=1,c=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | c<=0 | lambda<=0,
           NaN,lambda*dexpo.weibull(lambda*x,
           alpha=c,theta=alpha)/
           (1-pexpo.weibull(lambda*x,
           alpha=c,theta=alpha)))
	return(ret)
}

#The functions {\sf dexpo.weibull} and {\sf pexpo.weibull} are  from the {\sf R} contributed package {\sf reliaR}.



#Exponentiated logistic hazard rate function

hel=function(x,alpha=1,beta=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | beta<=0, NaN,
           dexpo.logistic(x,
           alpha=alpha,beta=beta)/
           (1-pexpo.logistic(x,
           alpha=alpha,beta=beta)))
	return(ret)
}

#The functions {\sf dexpo.logistic} and {\sf pexpo.logistic} are  from the {\sf R} contributed package {\sf reliaR}.


#BurrX hazard rate function

hburrx=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0, NaN,
           dburrX(x,alpha=alpha,
           lambda=lambda)/
           (1-pburrX(x,alpha=alpha,
           lambda=lambda)))
	return(ret)
}

#The functions {\sf dburrX} and {\sf pburrX} are  from the {\sf R} contributed package {\sf reliaR}.



#Gumbel hazard rate function

hgumbel=function(x,mu=1,sigma=1)
{
    ret	= ifelse(sigma<=0, NaN,
           dgumbel(x,mu=mu,sigma=sigma)/
           (1-pgumbel(x,mu=mu,sigma=sigma)))
	return(ret)
}

#The functions {\sf dgumbel} and {\sf pgumbel} are  from the {\sf R} contributed package {\sf reliaR}.


#Exponential extension hazard rate function

hexpext=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0,
           NaN, alpha*lambda*
           (1+lambda*x)**(alpha-1))
	return(ret)
}




#Chen hazard rate function

hchen=function(x,beta=1,lambda=1)
{
    ret	= ifelse(x<=0 | beta<=0 | lambda<=0,
           NaN, beta*lambda*
           x**(beta-1)*exp(x**beta))
	return(ret)
}





#Loggamma hazard rate function

hlgamma=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=1 | alpha<=0 | lambda<=0, NaN,
           dlgamma(x,shapelog=alpha,
           ratelog=lambda)/
           (1-plgamma(x,shapelog=alpha,
           ratelog=lambda)))
	return(ret)
}

#The functions {\sf dlgamma} and {\sf plgamma} are  from the {\sf R} contributed package {\sf actuar}.


#Marshall-Olkin exponential hazard rate function

hmoe=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0,
           NaN, lambda*
           exp(lambda*x)/
           (exp(lambda*x)-1+alpha))
	return(ret)
}




#Marshall-Olkin Weibull hazard rate function

hmow=function(x,alpha=1,beta=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | beta<=0 | lambda<=0,
           NaN, beta*
           lambda**beta*x**{beta-1}*
           exp((lambda*x)**beta)/
           (exp((lambda*x)**beta)-1+alpha))
	return(ret)
}




#Logistic exponential hazard rate function

hle=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0, NaN,
           dlogis.exp(x,
           alpha=alpha,lambda=lambda)/
           (1-plogis.exp(x,
           alpha=alpha,lambda=lambda)))
	return(ret)
}

#The functions {\sf dlogis.exp} and {\sf plogis.exp} are  from the {\sf R} contributed package {\sf reliaR}.



#Logistic Rayleigh hazard rate function

hlr=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0, NaN,
           dlogis.rayleigh(x,
           alpha=alpha,lambda=lambda)/
           (1-plogis.rayleigh(x,
           alpha=alpha,lambda=lambda)))
	return(ret)
}

#The functions {\sf dlogis.rayleigh} and {\sf plogis.rayleigh} are  from the {\sf R} contributed package {\sf reliaR}.



#Loglog hazard rate function

hll=function(x,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | lambda<=0,
           NaN,alpha*log(lambda)*
           x**(alpha-1)*
           lambda**(x**alpha))
	return(ret)
}




#Generalized power Weibull hazard rate function

hgpw=function(x,alpha=1,theta=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | theta<=0,
           NaN,alpha*theta*x**{alpha-1}*
           (1+x**alpha)**{theta-1})
	return(ret)
}




#Flexible Weibull hazard rate function

hfw=function(x,alpha=1,beta=1)
{
    ret	= ifelse(x<=0 | alpha<=0 | beta<=0,
           NaN,(alpha+beta/(x*x))*
           exp(alpha*x-beta/x))
	return(ret)
}




#Generalized $F$ hazard rate function

hgenF=function(x,beta=0,sigma=1,m1=1,m2=1)
{
    ret	= ifelse(x<=0 | sigma<=0 | m1<=0 | m2<=0,
           NaN,exp(-beta*m1/sigma)*(m1/m2)**m1*x**(m1/sigma-1)/
           (sigma*beta(m1,m2)*(1+(m1/m2)*(exp(-beta)*x)**(1/sigma))**(m1+m2)*
           (1-pbeta(m1*(exp(-beta)*x)**(1/sigma)/(m2+m1*(exp(-beta)*x)**(1/sigma)),shape1=m1,shape2=m2))))
	return(ret)
}

#The function {\sf pbeta} is  from the {\sf R} base package.



#Gompertz integrated hazard rate function

igompertz=function(x,t=1,alpha=1,beta=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) | alpha<=0 | beta<=0,
           NaN, (1/beta)*
           (exp(beta*t)-1)*
           alpha*exp(beta*x))
	return(ret)
}



#Makeham integrated hazard rate function

imakeham=function(x,t=1,alpha=1,beta=1,epsilon=1)
{
	ret	=ifelse(x<=0 | t<=0 | length(x)!=length(t) | alpha<=0 | beta<=0 | epsilon<=0,
          NaN,t*epsilon+(1/beta)*
          (exp(beta*t)-1)*
          alpha*exp(beta*x))
	return(ret)
}



#Perks integrated hazard rate function

iperks=function(x,t=1,alpha=1,beta=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) | alpha<=0 | beta<=0,
           NaN, (1/beta)*log(
           (1+alpha*exp(beta*(x+t))) /
           (1+alpha*exp(beta*x))))
	return(ret)
}



#Beard integrated hazard rate function

ibeard=function(x,t=1,alpha=1,beta=1,rho=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) | alpha<=0 | beta<=0 | rho<=0,
           NaN, (1/beta)*(1/rho)*
           log((1+alpha*rho*exp(beta*(x+t))) /
           (1+alpha*rho*exp(beta*x))))
	return(ret)
}



#Makeham-Perks integrated hazard rate function

imakehamperks=function(x,t=1,alpha=1,beta=1,epsilon=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) | alpha<=0 | beta<=0 | epsilon<=0,
           NaN, t*epsilon+(1/beta)*
           (1-epsilon)*
           log((1+alpha*exp(beta*(x+t))) /
           (1+alpha*exp(beta*x))))
	return(ret)
}



#Makeham-Beard integrated hazard rate function

imakehambeard=function(x,t=1,alpha=1,beta=1,rho=1,epsilon=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) | alpha<=0 | beta<=0 | rho<=0 | epsilon<=0,
           NaN, t*epsilon+(1/beta)*
           (1/rho-epsilon)*
           log((1+alpha*rho*exp(beta*(x+t))) /
           (1+alpha*rho*exp(beta*x))))
	return(ret)
}



#Exponential integrated hazard rate function

iexponential=function(x,t=1,alpha=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) | alpha<=0,
           NaN,t*alpha)
	return(ret)
}



#Pareto integrated hazard rate function

ipareto=function(x,t=1,alpha=1,m=1)
{
	ret	= ifelse(x<=m | t<=0 | m<=0 | length(x)!=length(t) | alpha<=0,
           NaN,alpha*log(1+t/x))
	return(ret)
}



#Weibull integrated hazard rate function

iweibull=function(x,t=1,alpha=1,sigma=1)
{
	ret=ifelse(x<=0 | t<=0 | length(x)!=length(t) | alpha<=0 | sigma<=0,
                   NaN,((x+t)**sigma-x**sigma)/alpha**sigma)
        return(ret)
}



#Logistic integrated hazard rate function

ilogistic=function(x,t=1,alpha=1,sigma=1)
{
	ret	= ifelse(length(x)!=length(t) | sigma<=0,NaN,
           log((1 + exp((x+t-alpha)/sigma)) /
           (1 + exp((x-alpha)/sigma))))
	return(ret)
}



#Log-logistic integrated hazard rate function

iloglogistic=function(x,t=1,alpha=1,sigma=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) | alpha<=0 | sigma<=0,
           NaN, log(
           (1+alpha*(x+t)**sigma) /
           (1+alpha*x**sigma)))
	return(ret)
}




#Normal integrated hazard rate function

inormal=function(x,t=1,alpha=1,sigma=1)
{
	ret	= ifelse(length(x)!=length(t) | sigma<=0,NaN,
           log((1-pnorm((x-alpha)/sigma)) /
           (1-pnorm((x+t-alpha)/sigma))))
	return(ret)
}

#The function {\sf pnorm} is  from the {\sf R} base package.


#Lognormal integrated hazard rate function

ilognormal=function(x,t=1,alpha=1,sigma=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) | sigma<=0,
           NaN, log(
           (1-plnorm((log(x)-alpha)/sigma)) /
           (1-plnorm((log(x+t)-alpha)/sigma))))
	return(ret)
}

#The function {\sf plnorm} is  from the {\sf R} base package.



#Inverse-Gaussian integrated hazard rate function

iinvgauss=function(x,t=1,alpha=1,sigma=1)
{
	ret=NaN
        if (x>0&t>0&length(x)==length(t)&alpha>0&sigma>0)
        {tt=pnorm(sqrt(sigma/x)*(x/alpha-1))+exp(2*sigma/alpha)*pnorm(-sqrt(sigma/x)*(x/alpha+1))
        tt2=pnorm(sqrt(sigma/(x+t))*((x+t)/alpha-1))+exp(2*sigma/alpha)*pnorm(-sqrt(sigma/(x+t))*((x+t)/alpha+1))
        ret=log((1-tt)/(1-tt2))}
	return(ret)
}

#The function {\sf pnorm} is from the {\sf R} base package.


#Gamma integrated hazard rate function

igamma=function(x,t=1,alpha=1,lambda=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) | alpha<=0 | lambda<=0,
           NaN, log((1-pgamma(x,shape=lambda,scale=alpha)) /
           (1-pgamma(x+t,shape=lambda,scale=alpha))))
	return(ret)
}

#The function {\sf pgamma} is  from the {\sf R} base package.


#Generalized-gamma integrated hazard rate function

igengamma=function(x,t=1,b=1,d=1,k=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) | b<=0 | d<=0 | k<=0,
           NaN, log((1-pgamma((x/b)**d,shape=k))/(1-pgamma(((x+t)/b)**d,shape=k))))
	return(ret)
}

#The function {\sf pgamma} is  from the {\sf R} base package.


#Linear integrated hazard rate function

ilinear=function(x,t=1,a=1,b=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) | a<0 | b<0 | a==0&b==0,
           NaN, (a+b*x)*t+b*t*t/2)
	return(ret)
}



#Uniform integrated hazard rate function

iuniform=function(x,t=1,a=0,b=1)
{
	ret	= ifelse(x<=a | x>=b | t<=0 | length(x)!=length(t) | b<=a | a<=0,
           NaN,log((b-x)/(b-x-t)))
	return(ret)
}



#Kumaraswamy integrated hazard rate function

ikum=function(x,t=1,a=0,b=1)
{
	ret	= ifelse(x<=0 | x>=1 | t<=0 | t>=1 |
           length(x)!=length(t) | a<=0 | b<=0,
           NaN,b*log((1-x**a)/(1-(x+t)**a)))
	return(ret)
}




#Gumbel II integrated hazard rate function

ifrechet=function(x,t=1,a=1,b=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) | a<=0 | b<=0,
           NaN, b*(x**(-a-1)-(x+t)**(-a-1)))
	return(ret)
}




#Hjorth integrated hazard rate function

ihjorth=function(x,t=1,delta=1,theta=1,beta=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t)
           | delta<=0 | theta<=0 | beta<=0,
           NaN, delta*x*t+
           delta*t*t/2+(theta/beta)*
           log(1+beta*t/(1+beta*x)))
	return(ret)
}



#Additive Weibull integrated hazard rate function

iaddweibull=function(x,t=1,a=1,b=1,c=1,d=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t)
           | a<=0 | b<=0 | c<=0 | d<=0,
           NaN, a**b*(x+t)**b+c**d*(x+t)**d)
	return(ret)
}




#Schabe integrated hazard rate function

ischabe=function(x,t=1,theta=1,gamma=0.5)
{
	ret	= ifelse(x>theta | t>theta | length(x)!=length(t)
           | theta<=0 | gamma<=0 | gamma>=1,
           NaN,log(1+t/(theta*gamma+x))-
           log(1-t/(theta-x)))
	return(ret)
}




#Lai integrated hazard rate function

ilai=function(x,t=1,lambda=1,beta=1,nu=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) |
           lambda<=0 | beta<0 | nu<0,NaN,lambda*(-x**beta*
           genhypergeo(U=beta,L=1+beta,z=nu*x)
           +(x+t)**beta*
           genhypergeo(U=beta,L=1+beta,z=nu*(x+t))
           -(nu/(beta+1))*x**(beta+1)*
           genhypergeo(U=1+beta,L=2+beta,z=nu*x)
           +(nu/(beta+1))*(x+t)**(beta+1)*
           genhypergeo(U=1+beta,L=2+beta,z=nu*(x+t))))
	return(ret)
}

#The function {\sf genhypergeo} is  from the {\sf R} contributed package {\sf hypergeo}.


#Xie integrated hazard rate function

ixie=function(x,t=1,lambda=1,alpha=1,beta=1)
{
	ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t)
           | lambda<=0 | alpha<=0 | beta<=0,
           NaN,lambda*alpha*
           (exp(((x+t)/alpha)**beta)-
           exp((x/alpha)**beta)))
	return(ret)
}




#$J$-shaped integrated hazard rate function

ijshape=function(x,t=1,b=1,nu=1)
{
	y1=1-x/b
    y2=1-(x+t)/b
    ret	= ifelse(y1<=0 | y1>=1 | t<=0 | t>=b
           | length(x)!=length(t) | b<=0 | nu<=0,
           NaN,log((1-(1-y1*y1)**nu)/
           (1-(1-y2*y2)**nu)))
	return(ret)
}




#Beta integrated hazard rate function

ibeta=function(x,t=1,a=1,b=1)
{
    ret	= ifelse(x<=0 | x>=1 | t<=0 | t>=1 |
           length(x)!=length(t) | a<=0 | b<=0,NaN,
           log((1-pbeta(x,shape1=a,shape2=b))/
           (1-pbeta(x+t,shape1=a,shape2=b))))
	return(ret)
}

#The function  {\sf pbeta} is  from the {\sf R} base package.

#Exponentiated exponential integrated hazard rate function

iee=function(x,t=1,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) | alpha<=0 | lambda<=0,
           NaN,-log((1-pgen.exp(x+t,
           alpha=alpha,lambda=lambda))/
           (1-pgen.exp(x,
           alpha=alpha,lambda=lambda))))
	return(ret)
}

#The function {\sf pgen.exp} is  from the {\sf R} contributed package {\sf reliaR}.


#Inverse generalized exponential integrated hazard rate function

iige=function(x,t=1,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) |
           alpha<=0 | lambda<=0,NaN,-log((1-pinv.genexp(x+t,
           alpha=alpha,lambda=lambda))/
           (1-pinv.genexp(x,
           alpha=alpha,lambda=lambda))))
	return(ret)
}

#The function {\sf pinv.genexp} is  from the {\sf R} contributed package {\sf reliaR}.




#Exponentiated Weibull integrated hazard rate function

iew=function(x,t=1,alpha=1,c=1,lambda=1)
{
    ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t)
           | alpha<=0 | c<=0 | lambda<=0,
           NaN,-log((1-pexpo.weibull(lambda*(x+t),
           alpha=c,theta=alpha))/
           (1-pexpo.weibull(lambda*x,
           alpha=c,theta=alpha))))
	return(ret)
}

#The function {\sf pexpo.weibull} is  from the {\sf R} contributed package {\sf reliaR}.



#Exponentiated logistic integrated hazard rate function

iel=function(x,t=1,alpha=1,beta=1)
{
    ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t)
           | alpha<=0 | beta<=0,
           NaN,log((1-pexpo.logistic(x,
           alpha=alpha,beta=beta))/
           (1-pexpo.logistic(x+t,
           alpha=alpha,beta=beta))))
	return(ret)
}

#The function {\sf pexpo.logistic} is  from the {\sf R} contributed package {\sf reliaR}.



#BurrX integrated hazard rate function

iburrx=function(x,t=1,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t)
           | alpha<=0 | lambda<=0, NaN,-log((1-pburrX(x+t,
           alpha=alpha,lambda=lambda))/
           1-pburrX(x,
           alpha=alpha,lambda=lambda)))
	return(ret)
}

#The function {\sf pburrX} is  from the {\sf R} contributed package {\sf reliaR}.


#Gumbel integrated hazard rate function

igumbel=function(x,t=1,mu=1,sigma=1)
{
    ret	= ifelse(length(x)!=length(t) | sigma<=0,NaN,
           -log((1-pgumbel(x+t,
           mu=mu,sigma=sigma))/
           (1-pgumbel(x,
           mu=mu,sigma=sigma))))
	return(ret)
}

#The function {\sf pgumbel} is  from the {\sf R} contributed package {\sf reliaR}.



#Exponential extension integrated hazard rate function

iexpext=function(x,t=1,alpha=1,lambda=1)
{
    ret	= ifelse(x<= 0 | t<=0 | length(x)!=length(t)
           | alpha<=0 | lambda<=0,NaN,
           (1+lambda*x+lambda*t)**alpha-
           (1+lambda*x)**alpha)
	return(ret)
}




#Chen integrated hazard rate function

ichen=function(x,t=1,beta=1,lambda=1)
{
    ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) |
           beta<=0 | lambda<=0,
           NaN, lambda*
           (exp((x+t)**beta)-exp(x**beta)))
	return(ret)
}


#Loggamma integrated hazard rate function

ilgamma=function(x,t=1,alpha=1,lambda=1)
{
    ret	= ifelse(x<=1 | t<=0 | length(x)!=length(t)
           | alpha<=0 | lambda<=0,
           -log((1-plgamma(x+t,shapelog=alpha,
           ratelog=lambda))/
           (1-plgamma(x,shapelog=alpha,
           ratelog=lambda))))
	return(ret)
}

#The function {\sf plgamma} is  from the {\sf R} contributed package {\sf actuar}.



#Marshall-Olkin exponential integrated hazard rate function

imoe=function(x,t=1,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t)
           | alpha<=0 | lambda<=0,
           NaN, log((exp(lambda*(x+t))-1+alpha)/
           (exp(lambda*x)-1+alpha)))
	return(ret)
}




#Marshall-Olkin Weibull integrated hazard rate function

imow=function(x,t=1,alpha=1,beta=1,lambda=1)
{
    ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t)
           | alpha<=0 | beta<=0 | lambda<=0,
           NaN, log(
           (exp((lambda*(x+t))**beta)-1+alpha)/
           (exp((lambda*x)**beta)-1+alpha)))
	return(ret)
}


#Logistic exponential integrated hazard rate function

ile=function(x,t=1,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t)
           | alpha<=0 | lambda<=0,
           NaN,log((1-plogis.exp(x,
           alpha=alpha,lambda=lambda))/
           (1-plogis.exp(x+t,
           alpha=alpha,lambda=lambda))))
	return(ret)
}

#The function {\sf plogis.exp} is  from the {\sf R} contributed package {\sf reliaR}.


#Logistic Rayleigh integrated hazard rate function

ilr=function(x,t=1,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t)
           | alpha<=0 | lambda<=0,
           NaN,log((1-plogis.rayleigh(x,
           alpha=alpha,lambda=lambda))/
           (1-plogis.rayleigh(x+t,
           alpha=alpha,lambda=lambda))))
	return(ret)
}

#The function {\sf plogis.rayleigh} is  from the {\sf R} contributed package {\sf reliaR}.


#Loglog integrated hazard rate function

ill=function(x,t=1,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t) |
           alpha<=0 | lambda<=0,
           NaN,lambda**((x+t)**alpha)-
           lambda**(x**alpha))
	return(ret)
}


#Generalized power Weibull integrated hazard rate function

igpw=function(x,t=1,alpha=1,theta=1)
{
    ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t)
           | alpha<=0 | theta<=0,
           NaN,(1+(x+t)**alpha)**theta-
           (1+x**alpha)**theta)
	return(ret)
}


#Flexible Weibull integrated hazard rate function

ifw=function(x,t=1,alpha=1,beta=1)
{
    ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t)
           | alpha<=0 | beta<=0,
           NaN,exp(alpha*(x+t)-beta/(x+t))-
           exp(alpha*x-beta/x))
	return(ret)
}



#Generalized $F$ integrated hazard rate function

igenF=function(x,t=1,beta=0,sigma=1,m1=1,m2=1)
{
    ret	= ifelse(x<=0 | t<=0 | length(x)!=length(t)
           | sigma<=0 | m1<=0 | m2<=0,
           NaN,log((1-pbeta(m1*(exp(-beta)*x)**(1/sigma)/(m2+m1*(exp(-beta)*x)**(1/sigma)),shape1=m1,shape2=m2))/
           (1-pbeta(m1*(exp(-beta)*(x+t))**(1/sigma)/(m2+m1*(exp(-beta)*(x+t))**(1/sigma)),shape1=m1,shape2=m2))))
	return(ret)
}

#The function {\sf pbeta} is  from the {\sf R} base package.




#Gompertz quantile function

qgompertz=function(x,u=0.5,alpha=1,beta=1)
{
	ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u) | alpha<=0 | beta<=0,
           NaN, (1/beta)*
           log(1-(beta/alpha)*exp(-beta*x)*log(u)))
	return(ret)
}



#Perks quantile function

qperks=function(x,u=0.5,alpha=1,beta=1)
{
	ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u) | alpha<=0 | beta<=0,
           NaN, (1/beta)*
           (log((u**(-beta))*
           (1+exp(alpha+beta*x))-1)-alpha)-x)
	return(ret)
}



#Beard quantile function

qbeard=function(x,u=0.5,alpha=1,beta=1,rho=1)
{
	ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u) | alpha<=0 | beta<=0 | rho<=0,
           NaN, (1/beta)*
           (log((u**(-beta*rho))*
           (1+alpha*rho*exp(beta*x))-1)-log(alpha)-log(rho))-x)
	return(ret)
}



#Exponential quantile function

qexponential=function(x,u=0.5,alpha=1)
{
	ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u) | alpha<=0,
           NaN,-log(u)/alpha)
	return(ret)
}



#Pareto quantile function

qpareto=function(x,u=0.5,alpha=1,m=1)
{
	ret	= ifelse(x<=m | u<=0 | u>=1 | m<=0 | length(x)!=length(u) | alpha<=0,
           NaN,x*exp(-log(u)/alpha)-x)
	return(ret)
}



#Weibull quantile function

qweibull=function(x,u=0.5,alpha=1,sigma=1)
{
	ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u) | alpha<=0 | sigma<=0,
           NaN, (x**sigma-alpha**sigma*log(u))**(1/sigma)-x)
	return(ret)
}



#Logistic quantile function

qlogistic=function(x,u=0.5,alpha=1,sigma=1)
{
	ret	= ifelse(u<=0 | u>=1 | length(x)!=length(u) | sigma<=0,
           NaN, sigma*
           (log(1+exp((x-alpha)/sigma)-u)-log(u))-x+alpha)
	return(ret)
}


#Log-logistic quantile function

qloglogis=function(x,u=0.5,alpha=1,sigma=1)
{
	ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u) | alpha<=0 | sigma<=0,
           NaN, exp((1/sigma)*
           (log((1/u)*(1+alpha*x**sigma)-1)-log(alpha))-x))
	return(ret)
}



#Normal quantile function

qnormal=function(x,u=0.5,alpha=1,sigma=1)
{
	ret	= ifelse(u<=0 | u>=1 | length(x)!=length(u) | sigma<=0,
           NaN, sigma*
           qnorm(1-u*(1-pnorm((x-alpha)/sigma)))+x*alpha)
	return(ret)
}

#The functions {\sf pnorm} and {\sf qnorm} are  from the {\sf R} base package.




#Lognormal quantile function

qlognormal=function(x,u=0.5,alpha=1,sigma=1)
{
	ret=NaN
    if (x==0 & u>0 & u<1 & length(x)==length(u) & sigma>0)
        ret=exp(sigma*qnorm(u)+alpha)
    if (x>0 & u>0 & u<1 & length(x)==length(u) & sigma>0)
       ret=exp(sigma*qnorm(1-u*
              (1-pnorm((log(x)-alpha)/sigma)))+alpha)
	return(ret)
}

#The functions {\sf pnorm} and {\sf qnorm} are  from the {\sf R} base package.


#Inverse-Gaussian quantile function

qinversegaussian=function(x,u=0.5,alpha=1,sigma=1)
{
	m=length(x)
        ret=rep(NaN,m)
        if (x>0&u>0&u<1&length(x)==length(u)&alpha>0&sigma>0)
        {ret=x
        for (i in 1:m)
        {ff=function (t)
        {tt=pnorm(sqrt(sigma/x[i])*(x[i]/alpha-1))+exp(2*sigma/alpha)*pnorm(-sqrt(sigma/x[i])*(x[i]/alpha+1))
        tt2=pnorm(sqrt(sigma/(x[i]+t))*((x[i]+t)/alpha-1))+exp(2*sigma/alpha)*pnorm(-sqrt(sigma/(x[i]+t))*((x[i]+t)/alpha+1))
        return(log((1-tt)/(1-tt2))+log(u[i]))}
        ret[i]=uniroot(ff,lower=0,upper=1000)$root}}
	return(ret)
}

#The function {\sf pnorm} is from the {\sf R} base package.



#Gamma quantile function

qgammad=function(x,u=0.5,alpha=1,lambda=1)
{
	ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u) | alpha<=0 | lambda<=0,
           NaN,qgamma(1-u*(1-pgamma(x,shape=lambda,scale=alpha)),shape=lambda,scale=alpha)-x)
	return(ret)
}

#The functions {\sf pgamma} and {\sf qgamma} are  from the {\sf R} base package.


#Generalized gamma quantile function

qgengammad=function(x,u=0.5,b=1,d=1,k=1)
{
	ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u) | b<=0 | d<=0 | k<=0,
           NaN,b*(qgamma(1-u*(1-pgamma((x/b)**d,shape=k)),shape=k))**(1/d)-x)
	return(ret)
}

#The functions {\sf pgamma} and {\sf qgamma} are  from the {\sf R} base package.

#Linear quantile function

qlinear=function(x,u=0.5,a=1,b=1)
{
	ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u) | a<0 | b<0 | a==0&b==0,
           NaN, (1/b)*(-(a+b*x)+
           sqrt(a*a+b*b*x*x+2*a*b*x+2*b*log(u))))
	return(ret)
}




#Uniform quantile function

quniform=function(x,u=0.5,a=1,b=2)
{
	ret	= ifelse(x<=a | x>=b | u<=0 | u>=1 | length(x)!=length(u) | b<=a | a<=0,
           NaN,(b-x)*(1-u))
	return(ret)
}



#Kumaraswamy quantile function

qkum=function(x,u=0.5,a=1,b=1)
{
	ret	= ifelse(x<=0 | x>=1 | u<=0 | u>=1 |
           length(x)!=length(u) | a<=0 | b<=0,
           NaN,(1-u**(-1/b)*
           (1-x**a))**(1/a)-x)
	return(ret)
}



#Gumbel II quantile function

qfrechet=function(x,u=0.5,a=1,b=1)
{
	ret	= ifelse(x<=0 | u<=0 | u>=1 |
           length(x)!=length(u) | a<=0 | b<=0,NaN,
           (x**{-a-1}+log(u)/b)**(-1/(a+1))-x)
	return(ret)
}




#Hjorth quantile function

qhjorth=function(x,u=0.5,delta=1,theta=1,beta=1)
{
	m=length(x)
    ret=rep(NaN,m)
    if (x>0&u>0&u<1&length(x)==length(u)&delta>0&theta>0&beta>0)
       {for (i in 1:m)
            {ff=function (t) {delta*x[i]*t+delta*t*t/2+
                     (theta/beta)*log(1+beta*t
                     /(1+beta*x[i]))+log(u[i])}
             ret[i]=uniroot(ff,lower=0,upper=100)$root}
       }
	return(ret)
}

#The function {\sf uniroot} is  from the {\sf R} base package.



#Additive Weibull quantile function

qaddweibull=function(x,u=0.5,a=1,b=1,c=1,d=1)
{
	m=length(x)
    ret=rep(NaN,m)
    if (x>0&u>0&u<1&length(x)==length(u)&a>0&b>0&c>0&d>0)
       {for (i in 1:m)
            {ff=function (t) {a**b*(x[i]+t)**b+
                     c**d*(x[i]+t)**d+log(u[i])}
             ret[i]=uniroot(ff,lower=0,upper=100)$root}
       }
	return(ret)
}

#The function {\sf uniroot} is  from the {\sf R} base package.


#Schabe quantile function

qschabe=function(x,u=0.5,theta=1,gamma=0.5)
{
	ret	= ifelse(x>theta | u<=0 | u>=1 | length(x)!=length(u)
           | theta<=0 | gamma<=0 | gamma>=1,NaN,
           theta-x-theta*(gamma+1)*
           (1-(theta*x)/(u*(theta-x))))
	return(ret)
}


#Lai quantile function

qlai=function(x,u=0.5,lambda=1,beta=1,nu=1)
{
	m=length(x)
    ret=rep(NaN,m)
    if (x>0&u>0&u<1&length(x)==length(u)&lambda>0&beta>=0&nu>=0)
       {for (i in 1:m)
            {ff=function (t) {lambda*(-x[i]**beta*
             genhypergeo(U=beta,L=1+beta,z=nu*x[i])
             +(x[i]+t)**beta*
             genhypergeo(U=beta,L=1+beta,z=nu*(x[i]+t))
             -(nu/(beta+1))*x[i]**(beta+1)*
             genhypergeo(U=1+beta,L=2+beta,z=nu*x[i])
             +(nu/(beta+1))*(x[i]+t)**(beta+1)*
             genhypergeo(U=1+beta,L=2+beta,z=nu*(x[i]+t)))}
             ret[i]=uniroot(ff,lower=0,upper=100)$root}
       }
	return(ret)
}

#The function {\sf uniroot} is  from the {\sf R} base package.
#The function {\sf genhypergeo} is  from the {\sf R} contributed package {\sf hypergeo}.



#Xie quantile function

qxie=function(x,u=0.5,lambda=1,alpha=1,beta=1)
{
	ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | lambda<=0 | alpha<=0 | beta<=0,NaN,
           alpha*(log(exp((x/alpha)**beta)-
           log(u)/(lambda*alpha)))**(1/beta)-x)
	return(ret)
}




#$J$-shaped quantile function

qjshape=function(x,u=0.5,b=1,nu=1)
{
	y=1-x/b
    ret	= ifelse(y<=0 | y>=1 | u<=0 | u>=1 | length(x)!=length(u)
           | b<=0 | nu<=0,NaN,b*(1-(1-(1-u*(1-(1-y*y)**nu)
           )**(1/nu))**(1/2))-x)
	return(ret)
}




#Beta quantile function

qbetad=function(x,u=0.5,a=1,b=1)
{
    ret	= ifelse(x<=0 | x>=1 | u<=0 | u>=1 | length(x)!=length(u)
           | a<=0 | b<=0,
           NaN,qbeta(1-u*(1-pbeta(x,
           shape1=a,shape2=b)),
           shape1=a,shape2=b)-x)
	return(ret)
}

#The functions  {\sf pbeta} and {\sf qbeta} are  from the {\sf R} base package.



#Exponentiated exponential quantile function

qee=function(x,u=0.5,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | alpha<=0 | lambda<=0,
           NaN,qgen.exp(1-u*(1-pgen.exp(x,
           alpha=alpha,lambda=lambda)),
           alpha=alpha,lambda=lambda)-x)
	return(ret)
}

#The functions {\sf pgen.exp} and {\sf qgen.exp} are  from the {\sf R} contributed package {\sf reliaR}.



#Inverse generalized  exponential quantile function

qige=function(x,u=0.5,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | alpha<=0 | lambda<=0,
           NaN,qinv.genexp(1-u*(1-pinv.genexp(x,
           alpha=alpha,lambda=lambda)),
           alpha=alpha,lambda=lambda)-x)
	return(ret)
}

#The functions {\sf pinv.genexp} and {\sf qinv.genexp} are  from the {\sf R} contributed package {\sf reliaR}.


#Exponentiated Weibull quantile function

qew=function(x,u=0.5,alpha=1,c=1,lambda=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | alpha<=0 | c<=0 | lambda<=0,NaN,
           (1/lambda)*qexpo.weibull(1-u*
           (1-pexpo.weibull(lambda*x,
           alpha=c,theta=alpha)),
           alpha=c,theta=alpha)-x)
	return(ret)
}

#The functions {\sf pexpo.weibull} and {\sf qpexpo.weibull} are  from the {\sf R} contributed package {\sf reliaR}.

#Exponentiated logistic quantile function

qel=function(x,u=0.5,alpha=1,beta=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | alpha<=0 | beta<=0,NaN,qexpo.logistic(1-u*
           (1-pexpo.logistic(x,
           alpha=alpha,beta=beta)),
           alpha=alpha,beta=beta)-x)
	return(ret)
}

#The functions {\sf pexpo.logistic} and {\sf qexpo.logistic} are  from the {\sf R} contributed package {\sf reliaR}.


#BurrX quantile function

qburrx=function(x,u=0.5,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | alpha<=0 | lambda<=0,NaN,qburrX(1-u*
           (1-pburrX(x,
           alpha=alpha,lambda=lambda)),
           alpha=alpha,lambda=lambda)-x)
	return(ret)
}

#The functions {\sf pburrX} and {\sf qburrX} are  from the {\sf R} contributed package {\sf reliaR}.


#Gumbel quantile function

qgumbeld=function(x,u=0.5,mu=1,sigma=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u) | sigma<=0,
           NaN,qgumbel(1-u*(1-pgumbel(x,
           mu=mu,sigma=sigma)),
           mu=mu,sigma=sigma)-x)
	return(ret)
}

#The functions {\sf pgumbel} and {\sf qgumbel} are  from the {\sf R} contributed package {\sf reliaR}.



#Exponential extension quantile function

qexpext=function(x,u=0.5,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | alpha<=0 | lambda<=0,NaN,
           (1/lambda)*(((1+lambda*x)**alpha-
           log(u))**(1/alpha)-1-lambda*x))
	return(ret)
}




#Chen quantile function

qchen=function(x,u=0.5,beta=1,lambda=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | beta<=0 | lambda<=0,NaN,
           (log(exp(x**beta)-
           log(u)/lambda))**(1/beta)-x)
	return(ret)
}




#Loggamma quantile function

qlgammad=function(x,u=0.5,alpha=1,lambda=1)
{
    ret	= ifelse(x<=1 | u<=0 | u>=1 | length(x)!=length(u)
           | alpha<=0 | lambda<=0,
           NaN,qlgamma(1-u*(1-plgamma(x,shapelog=alpha,
           ratelog=lambda)),shapelog=alpha,
           ratelog=lambda)-x)
	return(ret)
}

#The functions {\sf plgamma} and {\sf qlgamma} are  from the {\sf R} contributed package {\sf actuar}.



#Marshall-Olkin exponential quantile function

qmoe=function(x,u=0.5,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | alpha<=0 | lambda<=0,NaN,
           (1/lambda)*log(1-alpha+(1/u)*
           (exp(lambda*x)-1+alpha))-x)
	return(ret)
}




#Marshall-Olkin Weibull quantile function

qmow=function(x,u=0.5,alpha=1,beta=1,lambda=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | alpha<=0 | beta<=0 | lambda<=0,NaN,
           (1/lambda)*(log(1-alpha+(1/u)*
           (exp((lambda*x)**beta)-
           1+alpha)))**(1/beta)-x)
	return(ret)
}




#Logistic exponential quantile function

qle=function(x,u=0.5,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | alpha<=0 | lambda<=0,
           NaN,qlogis.exp(1-u*
           (1-plogis.exp(x,
           alpha=alpha,lambda=lambda)),
           alpha=alpha,lambda=lambda)-x)
	return(ret)
}

#The functions {\sf plogis.exp} and {\sf qlogis.exp} are  from the {\sf R} contributed package {\sf reliaR}.




#Logistic Rayleigh quantile function

qlr=function(x,u=0.5,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | alpha<=0 | lambda<=0,
           NaN,qlogis.rayleigh(1-u*
           (1-plogis.rayleigh(x,
           alpha=alpha,lambda=lambda)),
           alpha=alpha,lambda=lambda)-x)
	return(ret)
}

#The functions {\sf plogis.rayleigh} and {\sf qlogis.rayleigh} are  from the {\sf R} contributed package {\sf reliaR}.



#Loglog quantile function

qll=function(x,u=0.5,alpha=1,lambda=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | alpha<=0 | lambda<=0,NaN,
           ((-log(u)+lambda**(x**alpha))/
           (log(lambda)))**(1/alpha)-x)
	return(ret)
}




#Generalized power Weibull quantile function

qgpw=function(x,u=0.5,alpha=1,theta=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | alpha<=0 | theta<=0,NaN,
           (((1+x**alpha)**theta-
           log(u))**(1/theta)-1)**
           (1/alpha)-x)
	return(ret)
}




#Flexible Weibull quantile function

qfw=function(x,u=0.5,alpha=1,beta=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | alpha<=0 | beta<=0,NaN,
           (1/(2*alpha))*(log(-log(u)+
           exp(alpha*x-beta/x))+
           sqrt((log(-log(u)+
           exp(alpha*x-beta/x)))**2+
           4*alpha*beta))-x)
	return(ret)
}



#Generalized $F$ quantile function

qgenF=function(x,u=0.5,beta=0,sigma=1,m1=1,m2=1)
{
    ret	= ifelse(x<=0 | u<=0 | u>=1 | length(x)!=length(u)
           | sigma<=0 | m1<=0 | m2<=0,
           NaN,exp(beta)*m1**(-sigma)*(1/(1-qbeta(1-(1-pbeta(m1*(exp(-beta)*x)**(1/sigma)/
           (m2+m1*(exp(-beta)*x)**(1/sigma)),
           shape1=m1,shape2=m2))*u,shape1=m1,shape2=m2))-m2)**sigma-x)
	return(ret)
}

#The functions {\sf pbeta} and {\sf qbeta} are  from the {\sf R} base package.
