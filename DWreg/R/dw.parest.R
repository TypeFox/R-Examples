dw.parest<-function(data,method="likelihood",method.opt="NR")
{
par.lik<-function(x,method="NR")
{
loglik.dw<-function(par,x)
	{
	q<-par[1]
	beta<-par[2]
	if(q<=0 | q>=1 | beta<=0)
		NA
	else
		{
		y<-as.numeric(lapply(x,ddw,q=q,beta=beta))
      	lk<-sum(log(y))
		return(lk)
		}
	}
par.in<-c(par.prop(x)$q,par.prop(x)$beta)
a<-maxLik(loglik.dw,start=par.in,x=x,method=method)
res<-a$estimate
return(list(q=res[1],beta=res[2],hess=a$hessian))
}
par.prop<-function(data)
{
a.prop<-sum(data==0)/length(data) 
idx1<-sum(data==1)/length(data) 
q.est<-1-a.prop
if(q.est==1)
	stop('There are no zeros in your sample', call. = FALSE)
if(q.est==0)
	stop('There are only zeros in your sample', call. = FALSE)
if(idx1!=0)
	beta.est<-log(log(q.est-idx1)/log(q.est))/log(2)

k<-max(data)-1
k<-min(10000,k)

if(k>1)
{
idx.i<-NULL
for(i in 1: k){
	a.prop<-a.prop+sum(data==i)/length(data) 
	idx.i[i]<-1/log(i+1)*log(log(1-a.prop)/log(q.est))
	}
beta.est<-1/(k)*sum(idx.i)
}
return(list(q=q.est,beta=beta.est))
}
if(method=="proportion")
	res<-par.prop(data)
if(method=="likelihood")
	res<-par.lik(data,method=method.opt)
return(res)
}