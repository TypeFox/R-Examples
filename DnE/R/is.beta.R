is.beta <-
function(x,m,a=10,sita1=NULL,sita2=NULL)
{
	re=1;
	for(i in 1:length(x))
		if(x[i]<=0||x[i]>=1)
			re=-1;
	p=rep(0,m);
	y=rep(0,m);
	q=0;
	if(re==-1)
	{
		return(data.frame("state"=-1,"pvalue"=0));
	}
	else
	{
		x1=mean(x);sum=0;
		for(i in 1:length(x))
		{
			sum=sum+x[i]^2;
		}
		x2=sum/length(x);
		if(is.null(sita1)&&is.null(sita2))
		{
			sita1=x1^2*(1-x1)/(x2-x1^2)-x1;
			sita2=sita1*(1-x1)/x1;
			df=m-3;
		}
		else if(is.null(sita1))
		{
			sita1=x1*(x1-x2)/(x2-x1^2);
			df=m-2;
		}
		else if(is.null(sita2))
		{
			sita2=sita1*(x1-1);
			df=m-2;
		}
		else
		{
			df=m-1;
		}
		if(sita1>0)
		{
		for(i in 1:m)
		{
			p[i]=pbeta(i/m,sita1,sita2)-pbeta((i-1)/m,sita1,sita2);
			if(p[i]==0)
			{
				break;
			}
			for(j in 1:length(x))
				if(x[j]>=(i-1)/m && x[j]<i/m)
					y[i]=y[i]+1;
			q=q+(y[i]-(length(x)*p[i]))^2/(length(x)*p[i]);
		}
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
}
