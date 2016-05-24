library(tmvtnorm)
library(rgl)

# simulate x1, x2, x3 from truncated multivariate normal distribution
sigma = matrix(c(1, 0, 0, 
                 0, 1, 0,
                 0, 0, 1), 3, 3)

# not truncated
X  = rmvnorm(n=2000, mean=c(0,0,0), sigma=sigma)
# truncated
X2 = rtmvnorm(n=2000, mean=c(0,0,0), sigma=sigma, lower=c(-Inf,-Inf,-Inf), upper=c(0,1,Inf))

# display as 3D scatterplot
open3d()
plot3d(X[,1], X[,2], X[,3],    col="black", size=2, xlab=expression(x[1]), ylab=expression(x[2]), zlab=expression(x[3]))
plot3d(X2[,1], X2[,2], X2[,3], col="red", size=2, add=TRUE)

