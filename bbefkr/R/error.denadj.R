error.denadj <- function(ban, badj, eps, res.data)
{
    res = as.numeric(res.data)
    data.num = length(res)
    std.res = sd(res)
    epsilon = (res-mean(res))/std.res
    eps.std = (eps-mean(res))/std.res
    tem = tem2 = vector(,length(epsilon))
    for(i in 1:(length(epsilon)))
    {
        tem[i]=(eps.std-epsilon[i])/(ban*(1+badj*abs(res[i])))
        tem2[i]=dnorm(tem[i])/(ban*(1+badj*abs(res[i])))
    }		
    hatf=(sum(tem2)/data.num)/std.res
    return(hatf)		
}
