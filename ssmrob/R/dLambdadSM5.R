dLambdadSM5 <-
function(x, beta)  # for variance estimation  
{
  (-(1-pnorm(t(x)%*%beta))*dnorm(t(x)%*%beta)*(t(x)%*%beta)+
    dnorm(t(x)%*%beta)^2)/(1-pnorm(t(x)%*%beta))/(1-pnorm(t(x)%*%beta))
}
