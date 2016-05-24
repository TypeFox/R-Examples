ci.RateR <-
function(x1, x2, e1, e2, cov.prob, midp=T)
         {if(x1+x2==0) 
            {lower=0; upper=Inf; est=NA; p=1}
          if(x1+x2>0)
            {fit=binom.confint(x1, x1+x2, cov.prob)
             lower=e1*(1-fit$upper)/fit$upper/e2
             upper=e1*(1-fit$lower)/fit$lower/e2
             est=e1*x2/x1/e2
             p=2*min(pbinom(x1, x1+x2, e1/(e1+e2))-dbinom(x1, x1+x2, e1/(e1+e2))*0.5, 1-pbinom(x1, x1+x2, e1/(e1+e2))+dbinom(x1, x1+x2, e1/
(e1+e2))*0.5)
             }
           return(list(est=est, lower=lower, upper=upper, p=p))
          }
