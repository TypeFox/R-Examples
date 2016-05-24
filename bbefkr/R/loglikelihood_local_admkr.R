loglikelihood_local_admkr <- function(h, resid)
{
    b = h[2]
    badj = h[3]
    epsilon = scale(resid)
    std = sd(resid)
    cont = (2.0*pi)^(-0.5)
    logf = vector(,length(resid))
    for(i in 1:length(resid))
    {
        temp = epsilon[i] - epsilon[-i]
        res = sum(cont*exp(-0.5*((temp/(b*(1+badj*abs(epsilon[-i]))))^2))/(b*(1+badj*abs(epsilon[-i]))))
        logf[i] = log(res/length(temp)/std)
    }
    return(sum(logf))			
}
