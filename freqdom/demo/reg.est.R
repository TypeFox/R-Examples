library("freqdom")

# Parameters
d = 15      # dimension
n = 1500    # number of observation
ntest = 100 # test observations (not used for training)

# Generate process Y_t = P0 X_t + P1 X_{t-1} + e
P = 1/((1:d)%*%t(1:d))
P = P / norm.spec(P)

X = rar(n + ntest,Psi=matrix(0,d,d),noise=rnorm)
e = rar(n + ntest,Psi=matrix(0,d,d),noise=rnorm)

Y = X %*% P  + e * 0.1

# decaying operator
P = sqrt(d:1 %*% t(d:1))
P = P + matrix(rnorm(d*d),d,d)*0.0
P = P / norm.spec(P) # normalisation |P| = 1

par(mfrow=c(1,1))
persp(1:d,1:d,P,theta=190,phi=30, zlab="P[i,j]", main="Operator to estimate")

# generate n + ntest independent random vectors X, the same for noise
X = rar(n + ntest,Psi=matrix(0,d,d),noise=rnorm)
e = rar(n + ntest,Psi=matrix(0,d,d),noise=rnorm)

# compute the response
Y = X %*% P  + e * 0.1

# Estimate using LS method and data-driven dimension choise
P.est = reg.est(X[1:n,],Y[1:n,])

par(mfrow=c(1,2))
persp(1:d,1:d,P,theta=190,phi=30, zlab="P[i,j]", main="Operator to estimate")
persp(1:d,1:d,P.est,theta=190,phi=30, zlab="P[i,j]", main="LS estimator")

# Predict
Y.est = X %*% P.est

# Compare estimations and predictions
MSE(P,P.est)
MSE(P,0)

# in-sample
MSE(Y[1:n,],Y.est[1:n,])
MSE(Y[1:n,],0)

# out-sample
MSE(Y[n + 1:ntest,],Y.est[n + 1:ntest,])
MSE(Y[n + 1:ntest,],0)

# use spectral estimator
A = speclagreg(X,Y)
P.est.spec = A$operators[1,,]

par(mfrow=c(1,2))
persp(1:d,1:d,P,theta=190,phi=30, zlab="P[i,j]", main="Operator to estimate")
persp(1:d,1:d,P.est.spec,theta=190,phi=30, zlab="P[i,j]", main="Frequency domain estimator")

# Predict
Y.est = X %*% P.est.spec

# Compare estimations and predictions
MSE(P,P.est.spec)
MSE(P,0)

# in-sample
MSE(Y[1:n,],Y.est[1:n,])
MSE(Y[1:n,],0)

# out-sample
MSE(Y[n + 1:ntest,],Y.est[n + 1:ntest,])
MSE(Y[n + 1:ntest,],0)
par(mfrow=c(1,1))
