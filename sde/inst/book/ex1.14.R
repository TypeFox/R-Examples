# ex1.14.R -- corrected version. See errata corrige to the first edition
set.seed(123)
par("mar"=c(3,2,1,1))
par(mfrow=c(2,1))
npaths <- 30
N <- 1000
sigma <- 0.5
nu <- -0.7
X <- sde.sim(drift=expression(0),sigma=expression(0.5), pred=F, N=N,M=npaths) 
Y <- X + nu*time(X)
girsanov <- exp(0.25 * (nu/sigma*X[N,] + 0.5*(nu/sigma)^2))
girsanov <- (girsanov - min(girsanov)) / diff(range(girsanov))
col.girsanov <- gray(1-girsanov)
matplot(time(X),Y,type="l",lty=1, col="black",xlab="t")
matplot(time(X),Y,type="l",lty=1,col=col.girsanov,xlab="t")
