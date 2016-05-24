
'chesson' <-
  function(alpha=c(1.1*1e-5, 1e-5), d=.1, years=10, N0=c(1e3,1e5),
                    w=c(.6, 1), env.var=1,
                    specialization=1, spread=0.67, type=1)
{
  if(spread>1 | spread<0) {
    stop("'spread' must be between zero and 1.")
  }
  if(specialization < 0 ) {
    stop("'specialization' must be non-negative.")
  }
### In terms of the Beta distribution:
  ## 'specialization' is a - 1,
  ## b is a function of the mode which is a function of 'spread'.
t <- 1:years
env.beta <- rbeta(years, 1/env.var, 1/env.var)

w.rare <- w[1]
w.comm <- w[2]
a.r <- specialization + 1
mode.r <- .5 + spread/2
b.r <- (a.r - 1)/mode.r - a.r + 2
  if(is.infinite(a.r)) stop("Niche not logical (a of Beta(x,a,b) is infinite)")
a.c <- b.r
b.c <- a.r
Bs <- matrix(NA, nrow=years, ncol=2)
a.r;b.r
  Bs[,1] <- dbeta(env.beta, a.r, b.r)
Bs[,2] <- dbeta(env.beta, a.c, b.c)
rhof <- function(x) {min(dbeta(x, a.r,b.r),dbeta(x,a.c,b.c)) }
rho <- integrate(Vectorize(rhof), 0, 1)
Es <- matrix(NA, nrow=years, ncol=2)
Es[,1] <- w.rare * Bs[,1] 
Es[,2] <- w.comm * Bs[,2] 
alpha <- alpha
d <- d
Ns <- matrix(NA, nrow=years+1, ncol=2)
Ns[1,] <- N0
Cs <- matrix(NA, nrow=years, ncol=2)
Rs <- matrix(NA, nrow=years, ncol=2)

for(i in 1:years) Ns[i+1,] <- { 
  juveniles <- sum(exp(Es[i,])*Ns[i,])
  if(type==1) 
    Cs[i,]  <-  alpha*juveniles   else {
    Cs[i,] <- rep( log(juveniles/sum(d*Ns[i,])), 2) }
  Rs[i,] <- exp(Es[i,]-Cs[i,])
  (1-d) * Ns[i,] + Rs[i,]*Ns[i,]
                              }
 return(list(time=c(0,t), Ns=Ns, Es=Es, Cs=Cs, Rs=Rs, Bs=Bs,
             env=env.beta-.5, overlap=rho[["value"]],
             params=c(a.r,b.r,a.c,b.c, spread)))
}
