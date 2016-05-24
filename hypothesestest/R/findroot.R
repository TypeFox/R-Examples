findroot <-
function(alpha=0.05,side="both",method="n",n,mu=0,sigma=1)
{
	if(method=="chi")
	{
    z=numeric()
	f=function(x)(x^(n/2-1))*exp(-x/2)/(gamma(n/2)*2^(n/2))
	for(i in 1:2)
	{
	    a=100
	    b=0
	    aa=integrate(f,0,a)
        bb=integrate(f,0,b)
	    if(i==1) alp=alpha/2
	    else alp=1-alpha/2
	    while(abs(aa$value-alp)!=0&abs(bb$value-alp)!=0)
        { 
		    c=(a+b)/2
		    aa=integrate(f,0,a)
		    bb=integrate(f,0,b)
		    if(abs(aa$value-alp)<1e-6) {z[i]=a;break}
		    if(abs(bb$value-alp)<1e-6) {z[i]=b;break}
		    if(integrate(f,0,c)$value<alp) b=c 
            else a=c
        }
  	}
	return(z)
	}
    if(alpha>0.5|alpha<0.001) 
    {
         print("Alpha must between 0.001 to 0.5")
         return(-1)
    }
    if(method!="n"&method!="t"&method!="dn"&method!="dt") print("Please input n t dn dt or chi")
    if(side=="both") alpha=alpha/2
    if(method=="n"|method=="dn"){f=function(x)(1/(sqrt(2*pi)*sigma))*exp(-(x-mu)^2/(2*sigma^2));a=3;b=0}
    if(method=="t"|method=="dt"){f=function(x)gamma((1+n)/2)/(sqrt(pi*n)*gamma(n/2)*(1+x^2/n)^((n+1)/2));a=10;b=0}
    aa=integrate(f,-Inf,a)
    bb=integrate(f,-Inf,b)
    while(abs(aa$value+alpha-1)!=0&abs(bb$value+alpha-1)!=0)
    {
	  c=(a+b)/2
         aa=integrate(f,-Inf,a)
      	  bb=integrate(f,-Inf,b)
   	  if(abs(aa$value+alpha-1)<1e-6) {z=a;break}
	  if(abs(bb$value+alpha-1)<1e-6) {z=b;break}
	  if(integrate(f,-Inf,c)$value<1-alpha) b=c
	  else a=c
    }
    if(side=="left") z=-z
    return(z)
}
