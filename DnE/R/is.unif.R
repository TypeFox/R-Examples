is.unif <-
function(x,m,a=10,sita1=NULL,sita2=NULL)
{
	p=rep(0,m);
	y=rep(0,m);
	q=0;
	x1=mean(x);sum=0;
	for(i in 1:length(x))
	{
		sum=sum+x[i]^2;
	}
	x2=sum/length(x);
	if (3*x2-3*x1^2>=0)
	{
		if(is.null(sita1)&&is.null(sita2))
		{
			sita1=x1-sqrt(3*x2-3*x1^2);
			sita2=x1+sqrt(3*x2-3*x1^2);
			df=m-2;
		}
		else if(is.null(sita1))
		{
			sita1=x1-sqrt(3*x2-3*x1^2);
			df=m-1;
		}
		else if(is.null(sita2))
		{
			sita2=x1+sqrt(3*x2-3*x1^2);
			df=m-1;
		}
		else 
		{
			df=m;
		}
		di=max(x)-min(x);
		for(i in 1:m)
		{
			p[i]=punif(min(x)+di*i/m,min=sita1,max=sita2)-punif(min(x)+di*(i-1)/m,min=sita1,max=sita2);
			if(p[i]==0)
			{
				break;
			}
			for(j in 1:length(x))
				if(x[j]>(min(x)+di*(i-1)/m) && x[j]<=(min(x)+di*i/m))
					y[i]=y[i]+1;
			q=q+(y[i]-(length(x)*p[i]))^2/(length(x)*p[i]);
		}
                             p[m+1]=punif(Inf,sita1,sita2);
                             q=q+length(x)*p[m+1];
		q0=qchisq(1-a,df);
                pvalue=pchisq(q,df);
		if(q<=q0)
		{
			return(data.frame("qchisq"=q,"pvalue"=pvalue));
		}
		else
		{
			return(data.frame("state"=-1,"pvalue"=1));
		}
	}
	else
	{
		return(data.frame("state"=-1,"pvalue"=1));
	}
}
