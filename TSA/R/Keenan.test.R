`Keenan.test` <-
function(x,order,...){
if(missing(order)) order=ar(x,...)$order 
m=order
myfun=function(y,x,m){
X=NULL
for (i in 1:m) X=cbind(X,zlag(x,i))
X=cbind(y,X)
X=na.omit(X)
lm.fit(y=X[,1],x=X[,-1,drop=FALSE])
}
x=as.vector(x)
n=length(x)

lm1=myfun(y=x,x=x,m=m)
fit1=lm1$fit^2
res1=lm1$residuals

lm2=myfun(y=c(rep(NA,m),fit1),x=x,m=m)
res2=lm2$residuals

lm3=lm(res1~res2-1)
test.stat=summary(lm3)$fstatistic[1]*(n-2*m-2)/(n-m-1)
names(test.stat)=NULL
return(list(test.stat=test.stat, p.value=pf(test.stat,
df1=1,df2=n-2*m-2,lower.tail=FALSE),order=order))
}

