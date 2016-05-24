bprobit.probs=function(X1,fit)
{
# bprobit.probs  Produces a simulated sample from the posterior 
#             distribution of an expected response for a linear regression model
#         X1 = design matrix of interest
#         fit = output of bayes.probit function

d=dim(X1)
n1=d[1]
md=dim(fit); m=md[1]
m1=array(0,c(m,n1))

for (j in 1:n1)
{
m1[,j]=pnorm(X1[j,]%*%t(fit))
}
return(m1)
}






