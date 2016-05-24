conint <-
function(TrnX=NULL,TrnY=NULL,m,n1,n2,s1,s2,side="both",alpha=0.95,method="n")
{
       alpha=1-alpha;
	if(method!="n"&method!="t"&method!="dn"&method!="dt"&method!="chi") print("Please input n t dn dt or chi")
	if(method=="chi")
	{
		z=findroot(alpha,side,method,n1-1)
		a=z[1]
		b=z[2]
		print(n1*s1^2/b)
		print(n1*s1^2/a)
    }
    else
	{if((method=="dn"|method=="dt")&is.null(TrnY)==F) 
	{
		n1=length(TrnX)
		n2=length(TrnY)
		sn=n1+n2-2
		s1=sd(TrnX)*sqrt((n1-1)/n1)
		s2=sd(TrnY)*sqrt((n2-1)/n2)
		sm=mean(TrnX)-mean(TrnY)
	}
	if(method=="dn"|method=="dt") std=s1^2/n1+s2^2/n2
	if(method=="dt") z=findroot(alpha,side,method,n1+n2-2)
	else z=findroot(alpha,side,method,n1)
	if(is.null(TrnX)==T)
	{
		sm=m
		sn=n1
		std=s1
       }
       else if(method=="n"|method=="t")
	{
	    sm=mean(TrnX)
	    sn=length(TrnX)
           std=sd(TrnX)*sqrt((sn-1)/sn)
       }
	if(method=="n"|method=="dn") k=z*std/sqrt(sn)
	else if(method=="t") k=z*std/sqrt(sn-1)
	else k=z*sqrt((n1*s1^2+n2*s2^2)*(1/n1+1/n2)/(n1+n2-2))
	if(side=="right")
	{
           a=-Inf
	    b=sm+k
	}
	else if(side=="left")
	{
           a=sm-k
	    b=Inf
	}
	else
	{
           a=sm-k
	    b=sm+k
       }
       print(a)
       print(b)
    }
}
