library(fda)

n = 500 # number of observations to generate
basis = create.bspline.basis(c(0,1),25)

# Simulate X_t - n brownian motions
# Set Y_t = w BB(X_t) - (1-w)BB(X_(t-1))
# where BB is a brownian bridge
d = 100
w = 0.6 # weight of lag 0 (where 1-w is the weight for lag 1)
X = rbm(n+1,d=d)
Y = (1-w)*rbb(BM=X[-(n+1),])
Y = Y - w*rbb(BM=X[-1,])
X = X[-1,]

# map to functional data
grid = 0:(d-1)/(d-1)
X = center.fd(Data2fd(grid,t(X),basis))
Y = center.fd(Data2fd(grid,t(Y),basis))

plot.alpha(t(X$coef),rem=1:15)



# divide into training and test set
ntest = floor(n * 0.1)
n = n - ntest

# estimate on the training set
A = speclagreg(t(X[1:n]$coef), t(Y[1:n]$coef), lags=-5:5, K = 25)
W = reglag.significance(t(X[1:n]$coef), t(Y[1:n]$coef), A, alpha = 0.05, plot=TRUE)

# check accuracy on the test set
Yest = A %c% t(X$coef)
print(paste("Relative error: ", MSE(t(Y$coef[,1:ntest + n]),Yest[1:ntest + n,]) / MSE(t(Y$coef[,1:ntest + n]),0)))
