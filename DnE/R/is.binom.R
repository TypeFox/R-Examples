is.binom <-
function(x,m,a=10,n=NULL,p0=NULL)
{
	re=1;
	for(i in 1:length(x))
		if(x[i]<0||round(x[i])!=x[i])
			re=-1;
	p=rep(0,m);
	y=rep(0,m);
	q=0;
	if(re==-1)
	{
		return(data.frame("state"=-1,"pvalue"=1));
	}
	else
	{
		x1=mean(x);sum=0;
		for(i in 1:length(x))
		{
			sum=sum+x[i]^2;
		}
		x2=sum/length(x);
		if(is.null(n)&&is.null(p0))
		{
			p0=1+(x1^2-x2)/x1;
			n=round(x1/p0);
			df=m-3;
		}
		else if(is.null(n))
		{
			n=round(x1/p0);
			df=m-2;
		}
		else if(is.null(p0))
		{
			p0=1+(x1^2-x2)/x1;
			df=m-2;
		}
		else
		{
			df=m-1;
		}
		if(p0>0&&p0<1)
		{
		for(i in 1:m)
		{
			p[i]=pbinom(round(n*i/m),n,p0)-pbinom(round(n*(i-1)/m),n,p0);
			if(p[i]==0)
			{
				break;
			}
			for(j in 1:length(x))
				if(x[j]>n*(i-1)/m && x[j]<=n*i/m)
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
