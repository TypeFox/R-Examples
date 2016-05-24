discrete.bayes.2=function(df,prior,y=NULL,...)
{
like=function(i,...)
if(is.matrix(y)==TRUE)
df(y[i,],param1,param2,...) else
df(y[i],param1,param2,...)
n.rows=dim(prior)[1]
n.cols=dim(prior)[2]
param1=as.numeric(dimnames(prior)[[1]])
param2=as.numeric(dimnames(prior)[[2]])
param1=outer(param1,rep(1,n.cols))
param2=outer(rep(1,n.rows),param2)
likelihood=1
if(length(y)>0)
{
n=ifelse(is.matrix(y)==FALSE,length(y),dim(y)[1])
for(j in 1:n)
likelihood=likelihood*like(j,...)
}
product=prior*likelihood
pred=sum(prior*likelihood)
prob=prior*likelihood/pred
obj=list(prob=prob,pred=pred)
class(obj)<-"bayes2"
obj
}