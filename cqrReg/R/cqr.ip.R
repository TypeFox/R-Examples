cqr.ip=function(X,y,tau)
{

weight=rep(1/length(tau),length(tau))
a=rq.fit.hogg(X,y,tau,weight,matrix(0,dim(X)[1],dim(X)[2]+length(tau)),rep(0,length(y)))
intercept=a$coe[1:length(tau)]
beta=a$coe[(length(tau)+1):length(a$coe)]
names(intercept)=tau
return(list(intercepts=intercept,beta=beta))

}
