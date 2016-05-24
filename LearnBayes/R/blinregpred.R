blinregpred=function(X1,theta.sample)
{
#blinregpred  Produces a simulated sample from the posterior predictive
#             distribution of a linear regression model
#         X1 = design matrix of interest
#         theta.sample = output of blinreg function

d=dim(X1)
n1=d[1]
m=length(theta.sample$sigma)
y1=array(0,c(m,n1))

for (j in 1:n1)
{
y1[,j]=t(X1[j,]%*%t(theta.sample$beta))+rnorm(m)*theta.sample$sigma
}
return(y1)
}






