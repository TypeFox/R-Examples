poisson.confint <-
function(x, prob)
                {alpha=1-prob
                 n=length(x)
                 lower=upper=rep(0, n)
                 for(i in 1:n)
                    {x0=x[i]
                     if(x0==0)
                        {f=function(lambda){ppois(x0, lambda)-dpois(x0, lambda)*0.5-alpha}
                         fit=uniroot(f, c(0, 1000))
                         upper[i]=fit$root  
                         }
                     if(x0>0)
                        {f1=function(lambda){ppois(x0, lambda)-dpois(x0, lambda)*0.5-alpha/2}
                         f2=function(lambda){1-ppois(x0, lambda)+dpois(x0, lambda)*0.5-alpha/2}
                         fit1=uniroot(f1, c(0, 1000))
                         fit2=uniroot(f2, c(0, 1000))
                         upper[i]=fit1$root
                         lower[i]=fit2$root
                         }
                      }
                 return(list(lower=lower, upper=upper))
                 }
