require(tmvtnorm)
library(utils)

# Example 1 from Horrace (2005)
x1<-seq(-2, 3, by=0.1)
x2<-seq(-2, 3, by=0.1)

density<-function(x)
{
  sigma=matrix(c(1, -0.5, -0.5, 1), 2, 2)
  z=dtmvnorm(x, mean=c(0,0), sigma=sigma, lower=c(-1,-1))
  z
}

fgrid <- function(x, y, f)
{
    z <- matrix(nrow=length(x), ncol=length(y))
    for(m in 1:length(x)){
        for(n in 1:length(y)){
            z[m,n] <- f(c(x[m], y[n]))
        }
    }
    z
}

# compute the density function
d=fgrid(x1, x2, density)

# plot the density function as Contourplot
contour(x1, x2, d, nlevels=5, main="Truncated Multivariate Normal Density", xlab=expression(x[1]), ylab=expression(x[2]))
abline(v=-1, lty=3, lwd=2)
abline(h=-1,   lty=3, lwd=2)

# Example 2:
X=rtmvnorm(n=100, mean=c(0,0), sigma=matrix(c(1, 0.8, 0.8, 2), 2, 2), lower=c(-Inf,-Inf), upper=c(0,0))
plot(X, xlim=c(-3,3), ylim=c(-3,3), main="Samples from Multivariate Normal Distribution", xlab=expression(x[1]), ylab=expression(x[2]))
abline(v=0, lty=2)
abline(h=0, lty=2)

# Example 3: Profiling of rejection sampling: 10000 samples ~ 0.8 second
Rprof("rtmvnorm.out")
X=rtmvnorm(n=10000, mean=c(0,0), sigma=matrix(c(1, 0.8, 0.8, 2), 2, 2), lower=c(-Inf,-Inf), upper=c(0,0))
Rprof(NULL)
summaryRprof("rtmvnorm.out")

# Example 4: Profiling of Gibbs sampling: 10000 samples ~ 0.8 second
Rprof("rtmvnorm.gibbs.out")
m     = 10
a     = rep(-1, m)
b     = rep(1,  m)

# Erwartungswert und Kovarianzmatrix erzeugen
mu          = rep(0, m)
sigma       = matrix(0.8, m, m)
diag(sigma) = rep(1, m)

# Akzeptanzrate ausrechnen
alpha       = pmvnorm(lower=a, upper=b, mean=mu, sigma=sigma)
alpha

X=rtmvnorm(n=10000, mean=mu, sigma=sigma, lower=a, upper=b, algorithm="gibbs")
Rprof(NULL)
summaryRprof("rtmvnorm.gibbs.out")

# Sampling from non-truncated normal distribution 10000 samples ~ 0.02 second
Rprof("rmvnorm.out")
X=rmvnorm(n=10000, mean=c(0,0), sigma=matrix(c(1, 0.8, 0.8, 2), 2, 2))
Rprof(NULL)
summaryRprof("rmvnorm.out")
