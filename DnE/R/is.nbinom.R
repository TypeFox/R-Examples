is.nbinom <-
function(x,m,a=10,p0=NULL,r0=NULL)
{
	re=1;
	for(i in 1:length(x))
		if(x[i]<=0||round(x[i])!=x[i])
			re=-1;
	p=rep(0,m+2);
	y=rep(0,m+2);
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
		if(is.null(p0)&&is.null(r0))
		{
			p0=x1/(x2-x1^2);
			r0=round(p0*x1/(1-p0));
			df=m-1;
		}
		else if(is.null(p0))
		{
			p0=1-x1^2/(x2-x1^2);
			df=m;
		}
		else if(is.null(r0))
		{
			r0=round(p0*x1/(1-p0));
			df=m;
		}
		else
		{
			df=m+1;
		}
		if (p0>0&&p0<1&&r0>0)
		{
		di=max(x)-min(x);
		for(i in 1:m)
		{
			p[i]=pnbinom(min(x)+round(di*i/m),r0,p0)-pnbinom(min(x)+round(di*(i-1)/m),r0,p0);
			if(p[i]==0)
			{
				break;
			}
			for(j in 1:length(x))
				if(x[j]>(min(x)+di*(i-1)/m) &&x[j]<=(min(x)+di*i/m))
					y[i]=y[i]+1;
			q=q+(y[i]-(length(x)*p[i]))^2/(length(x)*p[i]);
		}
		p[m+1]=pnbinom(Inf,r0,p0)-pnbinom(max(x),r0,p0);
		p[m+2]=pnbinom(min(x),r0,p0);
		y[m+2]=length(which(x==min(x)));
		if(p[m+1]!=0)
			q=q+(y[m+1]-(length(x)*p[m+1]))^2/(length(x)*p[m+1]);
		if(p[m+2]!=0)
			q=q+(y[m+2]-(length(x)*p[m+2]))^2/(length(x)*p[m+2]);
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
