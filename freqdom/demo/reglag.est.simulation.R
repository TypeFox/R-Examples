# Parameters
d = 9      # dimension
n = 1000    # number of observation
ntest = 100 # test observations (not used for training)

# Generate AR(2) process
P = array(0,c(2,d,d))

P[1,,] = 1/((1:d)%*%t(1:d))
P[1,,] = 0.6 * P[1,,]/norm.spec(P[1,,])
P[2,,] = 1/((d:1)%*%t(1:d))
P[2,,] = 0.3 * P[2,,] / norm.spec(P[2,,])

A = timedom(P,0:1)

# generate n + ntest independent random vectors X, the same for noise
X = rar(n + ntest,Psi=matrix(0,d,d),noise=rnorm)

# compute the response
Y = linproc(X,A,noise = function(n){ rnorm(n) * 0.1 })

# Estimate using LS method and data-driven dimension choise
P.est = reg.est(X[1:n,],Y[1:n,])

# Predict
A.est = timedom(array(P.est,c(1,d,d)),0)
Y.est = linproc(X,A.est,noise = function(n){ rnorm(n) * 0 })

# in-sample
MSE(Y[1:n,],Y.est[1:n,])
MSE(Y[1:n,],0)

# out-sample
MSE(Y[n + 1:ntest,],Y.est[n + 1:ntest,])
MSE(Y[n + 1:ntest,],0)

# use spectral estimator
A.est = speclagreg(X,Y,lags=0:1)

par(mfrow=c(2,3))
persp(1:d,1:d,P[1,,],theta=190,phi=30, zlab="P[i,j]", main="A[0]")
persp(1:d,1:d,P.est,theta=190,phi=30, zlab="P[i,j]", main="LS A[0]")
persp(1:d,1:d,A.est$operators[1,,],theta=190,phi=30, zlab="P[i,j]", main="Spectral est A[0]")
persp(1:d,1:d,P[2,,],theta=190,phi=30, zlab="P[i,j]", main="A[1]")
D = matrix(0,d,d)
D[1,1]= 0.0001
persp(1:d,1:d,D,theta=190,phi=30, zlab="P[i,j]", main="LS A[1] = n/a (not yet)", zlim = c(0,1))
persp(1:d,1:d,A.est$operators[2,,],theta=190,phi=30, zlab="P[i,j]", main="Spectral est A[1]")

# Predict
Y.est = linproc(X,A.est,noise = function(n){ rnorm(n) * 0 })

# in-sample
MSE(Y[1:n,],Y.est[1:n,])
MSE(Y[1:n,],0)

# out-sample
MSE(Y[n + 1:ntest,],Y.est[n + 1:ntest,])
MSE(Y[n + 1:ntest,],0)

