binom.confint <-
function(x, n, alpha, method="exact")
      {if(x>0 && x<n)
         {
          error=(1-alpha)/2
          eq1=function(p){pbinom(x, n, p)-dbinom(x, n, p)*0.5-error}
          fit1=uniroot(eq1, c(0, 1))
          eq2=function(p){1-pbinom(x, n, p)+dbinom(x, n, p)*0.5-error}
          fit2=uniroot(eq2, c(0, 1))
          lower=fit2$root
          upper=fit1$root
          }
       if(x==0)
         {
          error=(1-alpha)
          eq1=function(p){pbinom(x, n, p)-dbinom(x, n, p)*0.5-error}
          fit1=uniroot(eq1, c(0, 1))
          lower=0
          upper=fit1$root

          }

       if(x==n)
         {error=(1-alpha)
          eq2=function(p){1-pbinom(x, n, p)+dbinom(x, n, p)*0.5-error}
          fit2=uniroot(eq2, c(0, 1))
          lower=fit2$root
          upper=1
          }

       return(list(lower=lower, upper=upper))
       }
