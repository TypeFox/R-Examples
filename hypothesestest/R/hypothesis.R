hypothesis <-
function (TrnX=NULL,TrnY=NULL,m,u0,n1,n2,s1=NULL,s2=NULL,sigma1=NULL,sigma2=NULL,alpha=0.05,method="n",H0="u=u0",p)
{
	if(H0=="u=u0") side="both"
	else if(H0=="u>u0"|H0=="u>=u0") side="left"
	else side="right"
	if(method!="chi1"&method!="chi2"&method!="chi3")
	{met="other"
       if(method=="n") 
       {
	     if(is.null(TrnX)==T)
	     {
		    sm=m
            if(is.null(sigma1)==F) 
            {
				std=sigma1
                t=findroot(alpha,side,method,n1)
                sn=n1
            }
            else 
	     {
		  method="t"
                std=s1
                t=findroot(alpha,side,method,n1-1)
                sn=n1-1
            }
	 }
	    else
		{
			sm=mean(TrnX)
			n1=length(TrnX)
			if(is.null(sigma1)==F)  
            {
				std=sigma1
                t=findroot(alpha,side,method,n1)
				sn=n1
			}
            else 
			{ 
				method="t"
				std=sd(TrnX)*sqrt((n1-1)/n1)
			       t=findroot(alpha,side,method,n1-1)
				sn=n1-1
			}

	    }
	    qq=(sm-u0)/(std/sqrt(sn))
	}
    if(method=="dn"&is.null(TrnY)==F)
	{
              n1=length(TrnX)
		n2=length(TrnY)
		s1=sd(TrnX)
		s2=sd(TrnY)
        sm=mean(TrnX)-mean(TrnY)
        if(is.null(sigma1)==T&is.null(sigma2)==T) 		
        {
			method="t"
			t=findroot(alpha,side,method,n1+n2-2)
			sn=n1+n2-2
            a=(n1*(s1)^2+n2*(s2)^2)/(n1+n2-2)
            qq=sm/(sqrt(a*(1/n1+1/n2)))
        }
        else
        {
			t=findroot(alpha,side,method,n1+n2)
			sn=n1+n2
            a=(sigma1)^2/n1+(sigma2)^2/n2
            qq=sm/(sqrt(a))
	    }

    }	
    if(method=="dn"&is.null(TrnY)==T)
	{
		sm=m
        if(is.null(sigma1)==T&is.null(sigma2)==T) 		
        {
			method="t"
			t=findroot(alpha,side,method,n1+n2-2)
			sn=n1+n2-2
            a=(n1*(s1)^2+n2*(s2)^2)/(n1+n2-2)
            qq=sm/(sqrt(a*(1/n1+1/n2)))
        }
        else 
        {
			t=findroot(alpha,side,method,n1+n2)
			sn=n1+n2
            a=(sigma1)^2/n1+(sigma2)^2/n2
            qq=sm/(sqrt(a))
	    }
	}
    if(H0=="u=u0")
	{ 
		if(abs(qq)>t) print("refuse H0")
        else print("we can not reject H0.")
    }
	else if(H0=="u<u0"|H0=="u<=u0")
    {
		if(qq>t) print("refuse H0")
        else print("we can not reject H0.")
    }
	else
    {
		if(qq<t) print("refuse H0")
        else print("we can not reject H0.")
    }
	}
    else
	{met="chi"
	if(method=="chi1")
	{
		
		if(is.null(TrnY)==F){n1=length(TrnX);p=TrnY/n1}
		qq=0
		for(i in 1:n1) qq=qq+(TrnX[i]-n1*p[i])^2/(n1*p[i])
		t=findroot(2*alpha,side,met,n1-1)[2]
		sn=n1-1
    }
    else if(method=="chi2")
	{
		n=numeric()
		n[1]=n1
		n[2]=n2
		ni=length(TrnX)
		X=data.frame(TrnX,TrnY)
		qq=0
		for(j in 1:2) for(i in 1:ni) qq=qq+(X[i,j]-n[j]*(X[i,1]+X[i,2])/(sum(n)))^2/(n[j]*(X[i,1]+X[i,2])/(sum(n)))
		t=findroot(2*alpha,side,met,ni-1)[2]
		sn=ni-1
	}
    else if(method=="chi3")
	{
		if(is.data.frame(TrnX)==F) {print("error,TrnX must be a data.frame");return(-1);}
		a=nrow(TrnX)
		b=ncol(TrnX)
		qq=0
		for(j in 1:b) for(i in 1:a) qq=qq+(X[i,j]-sum(X[i,])*sum(X[,j])/n1)^2/(sum(X[i,])*sum(X[,j])/n1)
		t=findroot(2*alpha,side,met,(a-1)*(b-1))[2]
		sn=(a-1)*(b-1)
	}
    if(qq>t) print("refuse H0")
	else print("we can not reject H0.")
    }
	print("t is");print(t)
	print("Q is");print(qq)
	if(method=="n"|method=="dn") f=function(x,sigma=1,mu=0)(1/(sqrt(2*pi)*sigma))*exp(-(x-mu)^2/(2*sigma^2))
    else if(method=="t") f=function(x,sigma=1,mu=0)gamma((1+sn)/2)/(sqrt(pi*sn)*gamma(sn/2)*(1+x^2/sn)^((sn+1)/2))
	else if(met=="chi") f=function(x,sigma=1,mu=0)(x^(sn/2-1))*exp(-x/2)/(gamma(sn/2)*2^(sn/2))
	if(side=="both"&met=="chi") pvalue=integrate(f,0,qq)
	else if(side=="both") pvalue=integrate(f,-abs(qq),abs(qq))
	else if(side=="left") pvalue=integrate(f,qq,Inf)
	else pvalue=integrate(f,-Inf,qq)
	print("p-value is");print(1-pvalue$value)
}
