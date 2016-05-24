library("freqdom")
library("MASS")

# Parameters
d = 15      # dimension
n = 200    # number of observation
ntest = 50 # test observations (not used for training)

# Set a operator on lag 0 and 1
P = array(0,c(2,d,d))
P[1,,] = 1/((1:d)%*%t(1:d))^2
P[1,,] = 0.3 * P[1,,]/norm.spec(P[1,,])
P[2,,] = 1/((d:1)%*%t(1:d))^2
P[2,,] = 0.6 * P[2,,] / norm.spec(P[2,,])
A = timedom(P,0:1)

# generate n + ntest independent random vectors X, the same for noise
#Psi = matrix(rnorm(d*d),d,d)
#Psi = 0.7*Psi / norm(Psi)
X = rar(n + ntest,Psi=matrix(0,d,d),noise=function(d){ 0.5*rmvnorm(1,rep(0,d)) })

# compute the response
Y = linproc(X,A,noise = function(n){ rnorm(n) * 0.2 })

A = speclagreg(X, Y, lags=-5:5)
R = reglag.boot(X, Y, A, rep=2, plot=TRUE)
A = timedom.trunc(A,R$suglags)

W = reglag.significance(X, Y, A, alpha = 0.05, plot=TRUE)
