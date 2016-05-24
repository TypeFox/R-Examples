is.gamma <-
function(x,m,a=10,a0=NULL,b0=NULL)
{
	re=1;
	for(i in 1:length(x))
		if(x[i]<0)
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
		x1=mean(x);
		x2=var(x)
		if(is.null(a0)&&is.null(b0))
		{
			b0=x1/x2;
			a0=x1*b0;
			df=m-1;
		}
		else if(is.null(a0))
		{
			a0=x1*b0;
			df=m;
		}
		else if(is.null(b0))
		{
			b0=x1/x2;
			df=m;
		}
		else
		{
			df=m+1;
		}
		if(a0>0&&b0>0)
		{
		di=max(x)-min(x);
		for(i in 1:m)
		{
			p[i]=pgamma(min(x)+di*i/m,a0,rate=b0)-pgamma(min(x)+di*(i-1)/m,a0,rate=b0);
			if(p[i]==0)
			{
				break;
			}
			for(j in 1:length(x))
			if(x[j]>(min(x)+di*(i-1)/m) && x[j]<=(min(x)+di*i/m))
					y[i]=y[i]+1;
			q=q+(y[i]-(length(x)*p[i]))^2/(length(x)*p[i]);
		}
		p[m+1]=pgamma(Inf,a0,rate=b0)-pgamma(max(x),a0,rate=b0);
		y[m+1]=0;
		p[m+2]=pgamma(min(x),a0,rate=b0);
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
