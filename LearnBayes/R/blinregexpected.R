blinregexpected=function(X1,theta.sample)
{
#blinregpred  Produces a simulated sample from the posterior 
#             distribution of an expected response for a linear regression model
#         X1 = design matrix of interest
#         theta.sample = output of blinreg function

d=dim(X1)
n1=d[1]
m=length(theta.sample$sigma)
m1=array(0,c(m,n1))

for (j in 1:n1)
{
m1[,j]=t(X1[j,]%*%t(theta.sample$beta))
}
return(m1)
}






