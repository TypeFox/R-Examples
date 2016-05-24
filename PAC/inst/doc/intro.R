## ---- echo=F-------------------------------------------------------------
set.seed(1)

## ------------------------------------------------------------------------
# Problem parameters
n = 5e3                       # number of observations
p = 1                         # number of dimensions
K = 3                         # number of clusters
w = rep(1,K)/K                # component weights
mu <- c(0,2,4)                # component means
sd <- rep(1,K)/K              # component standard deviations

# Generate K mixture of gaussian 
g <- sample(1:K,prob=w,size=n,replace=TRUE)   # ground truth for clustering
X <- as.matrix(rnorm(n=n,mean=mu[g],sd=sd[g]))

## ------------------------------------------------------------------------
library(PAC)
y <- PAC(X, K)

## ------------------------------------------------------------------------
print(fmeasure(g,y))

## ------------------------------------------------------------------------
y.deep <- PAC(X, K, maxlevel = 100)
print(fmeasure(g,y.deep))

## ------------------------------------------------------------------------
y.ll <- PAC(X, K, method = 'll')
print(fmeasure(g,y.ll))

