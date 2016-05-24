#######################################################
# Section 3.6 A Bayesian Test of the Fairness of a Coin
#######################################################

library(LearnBayes)

 pbinom(5, 20, 0.5)

 n = 20
 y = 5
 a = 10
 p = 0.5
 m1 = dbinom(y, n, p) * dbeta(p, a, a)/dbeta(p, a + y, a + n - 
     y)
 lambda = dbinom(y, n, p)/(dbinom(y, n, p) + m1)
 lambda

 pbetat(p,.5,c(a,a),c(y,n-y))

prob.fair=function(log.a)
{
 a = exp(log.a)
 m2 = dbinom(y, n, p) * dbeta(p, a, a)/
             dbeta(p, a + y, a + n - y)
 dbinom(y, n, p)/(dbinom(y, n, p) + m2)
}

n = 20; y = 5; p = 0.5
curve(prob.fair(x), from=-4, to=5, xlab="log a", 
  ylab="Prob(coin is fair)", lwd=2)

S=readline(prompt="Type  <Return>   to continue : ")

 n=20
 y=5
 a=10
 p=.5
 m2=0
 for (k in 0:y)
   m2=m2+dbinom(k,n,p)*dbeta(p,a,a)/dbeta(p,a+k,a+n-k)
 lambda=pbinom(y,n,p)/(pbinom(y,n,p)+m2)
 lambda
