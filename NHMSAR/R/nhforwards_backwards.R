nhforwards_backwards <-
function(prior, transition, obslik, filter_only=0) {


T = dim(obslik)[2]
Q1 = length(prior)

scale = matrix(1,1,T)
loglik = 0
alpha = matrix(0,Q1,T)
gamma = matrix(0,Q1,T)
xi = array(0,c(Q1,Q1,T-1))

t = 1
alpha[,1] = prior * obslik[,t]

scale[t] = sum(alpha[,t]) 
alpha[,t]  = normalise(alpha[,t])

for (t in 2:T) {
   tmp = (t(transition[,,t]) %*% alpha[,t-1]) * obslik[,t]
   scale[t] = sum(tmp)
   alpha[,t] = normalise(tmp)
   }
loglik = sum(log(scale))
beta = matrix(0,Q1,T)
gamma = matrix(0,Q1,T)
beta[,T] = matrix(1,Q1,1)
gamma[,T] = normalise(alpha[,T] * beta[,T])

for (t in seq(T-1,1,-1)) {
  b = beta[,t+1] * obslik[,t+1] 
  beta[,t] = normalise((transition[,,t+1] %*% b))
  gamma[,t] = normalise(alpha[,t] * beta[,t])
  xi[,,t] = normalise((transition[,,t+1] * (alpha[,t] %*% t(b))))
  }

FB = NULL
FB$gamma = gamma
FB$xi = xi
FB$loglik = loglik
FB$M = Q1
FB$alpha=alpha
FB$beta=beta
return(FB)

}
