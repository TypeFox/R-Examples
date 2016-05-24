
### test numeric methods, in particular handling of unequal
### function lengths
library(maxLik)

f <- function(x) {
   if(x[1] <= 0)
       return(NA)
                           # support of x[1] is (0, Inf)
   return(c(log(x[1]),x[2]))
}

ng <- numericGradient(f, c(0.01,1), eps=0.1)

nh <- try(numericHessian(f, t0=c(0.01,1), eps=0.1))
