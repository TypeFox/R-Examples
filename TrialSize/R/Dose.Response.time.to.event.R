Dose.Response.time.to.event <-
function(alpha, beta,T0,T,Ti,ci,fi){
lambdai=log(2)/Ti
lambda_mean=mean(lambdai)
epsilon=sum(ci*lambdai)
sigma<-function(lambdai)
{
sigma2=sqrt(lambdai^2*(1+(exp(-lambdai*T)*(1-exp(lambdai*T0)))/(T0*lambdai))^(-1))
}
sigma2=sigma(lambdai)
sigma0=sigma(lambda_mean)
n=(qnorm(1-alpha)*sigma0*sqrt(sum(ci^2/fi))+qnorm(1-beta)*sqrt(sum(ci^2*sigma2^2/fi)))^2/(epsilon^2)
}
