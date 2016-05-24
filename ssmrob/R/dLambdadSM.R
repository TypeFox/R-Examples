dLambdadSM <-
function(x, beta)  # for variance estimation	
{
  (-pnorm(t(x)%*%beta)*dnorm(t(x)%*%beta)*(t(x)%*%beta)-
    dnorm(t(x)%*%beta)^2)/pnorm(t(x)%*%beta)/pnorm(t(x)%*%beta)
}
