is.norm <-
function(x,m,a=10)
{
	x=scale(x);
	p=rep(0,m+2);
	y=rep(0,m+2);
	q=0;
		df=m+1;
		di=max(x)-min(x);
		for(i in 1:m)
		{
			p[i]=pnorm(di*i/m)-pnorm(di*(i-1)/m);
			if(p[i]==0)
			{
				break;
			}
			for(j in 1:length(x))
				if(x[j]>di*(i-1)/m && x[j]<=di*i/m)
					y[i]=y[i]+1;
			q=q+(y[i]-(length(x)*p[i]))^2/(length(x)*p[i]);
		}
		p[m+1]=pnorm(Inf)-pnorm(max(x));
		p[m+2]=pnorm(min(x));
		y[m+2]=length(which(x==min(x)));
		q=q+length(x)*p[m+1];
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
