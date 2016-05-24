## Comparing particle filtering with the Kalman filter
# Generate data from a local level model
N = 1000; W = 0.1^2; V = 1; m0 = 0; C0 = 1
true.x = rep(NA,N); true.x[1] = rnorm(1,m0,sqrt(C0))
for (i in 2:N) true.x[i] = rnorm(1,true.x[i-1],sqrt(W)) # Evolve x
y = rnorm(N,true.x,sqrt(V))                             # Noisy data

# For particle filter
J = 1e2
x  = matrix(NA,N,J); x[1,]  = rnorm(J,m0,C0) # Sample from the prior for x
ws = matrix(NA,N,J); ws[1,] = renormalize(dnorm(y[1],x[1,],sqrt(V),log=TRUE), log=T)

# For Kalman filter
m = rep(NA,N); m[1] = m0 # Kalman filter expectation
M = rep(NA,N); M[1] = C0 # Kalman filter variance

# Run both particle filter and Kalman filter
pb = txtProgressBar(0,N, style=3)
for (i in 2:N) {
  setTxtProgressBar(pb, i)

  # Particle filter
  component   = resample(ws[i-1,],J,"stratified","ess",0.8, log=F)
  x[i,]       = rnorm(J,x[i-1,component$indices],sqrt(W))
  log.weights = log(component$weights)+dnorm(y[i],x[i,],sqrt(V),log=TRUE)
  ws[i,]      = renormalize(log.weights,log=TRUE)

  # Kalman filter
  K    = (M[i-1]+W)/(M[i-1]+W+V) # Adaptive coefficient
  m[i] = K*y[i]+(1-K)*m[i-1]
  M[i] = K*V
}

pf.m = apply(x*ws,1,sum)
require(Hmisc)
pf.lb = pf.ub = numeric(N)
for (i in 1:N) 
{
  pf.lb[i] = wtd.quantile(x[i,], ws[i,], normwt=TRUE, probs=c(.025))
  pf.ub[i] = wtd.quantile(x[i,], ws[i,], normwt=TRUE, probs=c(.975))
}

plot(m,type='l',xlab='t',ylab='x', lty=2, ylim=range(pf.ub,pf.lb), 
     main="Means and 95% Credible Intervals for Kalman and Particle filter")
lines(m+qnorm(.975)*sqrt(M))
lines(m-qnorm(.975)*sqrt(M))
lines(pf.m,col='red', lty=2)
lines(pf.ub, col='red')
lines(pf.lb, col='red')
legend("bottomleft",inset=0.01,c("Kalman filter","Particle filter"),
       col=c("black","red"),lty=rep(1,2),bg="white")

