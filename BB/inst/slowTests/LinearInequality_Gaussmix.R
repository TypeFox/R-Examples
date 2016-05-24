
require(BB)

##########################
# Gaussian mixture density
dgaussmix <- function (p) {
prop <- p[1:nmix]
mu <- p[(nmix+1):(2*nmix)]
sigma <- p[(2*nmix+1)]
sapply(y, function(y)sum(prop*dnorm(y,mean=mu,sd=sqrt(sigma))))	
}

# generating random numbers from a Gaussian mixture
rgaussmix <- function (n, prop, mu, sigma) {
nmix <- length(mu)
imix <- sample(1:nmix, size=n, prob=prop, rep=TRUE)
y <- rnorm(n, mean = mu[imix], sd = sqrt(sigma))
return(y)
}

# Gaussian mixture  minus log-likelihood
gaussmix.mloglik <- function(p){
- sum(log(dgaussmix(p)))
}

# Gradient of Gaussian mixture log-likelihood
gaussmix.grad <- function(p){
g <- rep(NA, length(p))
f <- dgaussmix(p)
pj <- p[1:nmix]
mu <- p[(nmix+1): (2*nmix)]
sigma <- p[2*nmix + 1]
phi <- outer(y, mu, function(y, mu) dnorm(y,mean=mu,sd=sqrt(sigma)))
g[1:nmix] <- - colSums(phi/f)  
phi2 <- outer(y, mu, function(y, mu) (y - mu)/sigma)
fimuj <- t(t(phi * phi2) * pj)
g[(nmix+1): (2*nmix)] <- - colSums(fimuj/f)  
phi3 <- outer(y, mu, function(y, mu) (y - mu)^2/sigma)
fisig <- apply(t(t(phi * ( 1 - phi3) ) * pj), 1, sum)
g[2*nmix+1] <- sum(fisig / f) / (2 * sigma)
g
}

heq <- function(x) {
x[1] + x[2] + x[3] + x[4] - 1
}

hin <- function(x) {
h <- rep(NA, 9)
h[1] <- x[1] 
h[2] <- x[2] 
h[3] <- x[3] 
h[4] <- x[4]
h[5] <- 1 - x[1] 
h[6] <- 1 - x[2] 
h[7] <- 1 - x[3] 
h[8] <- 1 - x[4]
h[9] <- x[9] 
h
}

Amat <- matrix(0, 10, 9)
Amat[1, 1:4] <- 1  # corresponds to equality
Amat[2,1] <- Amat[3,2] <- Amat[4,3] <- Amat[5,4] <- Amat[10,9] <- 1
Amat[6, 1] <- Amat[7, 2] <- Amat[8, 3] <- Amat[9, 4] <- -1

b <- c(1,0,0,0,0,-1,-1,-1,-1, 0)
meq <- 1

# A data realization
p <- c(0.2,0.4,0.2,0.2)
nmix <- length(p)
mu <- c(0,3,7,11)
sigma <- 2
npts <- 500
set.seed(12345)
y <- rgaussmix(npts, p, mu, sigma)


ymean <- mean(y)
ysd <- sd(y)
p0 <- rep(1/nmix, nmix)
ymean0 <- ymean + ysd * runif(nmix, -1.2, 1.2)
ysd0 <- ysd
par0 <- c(p0,ymean0, ysd0^2)


# # The inequalities are defined such that:  Amat %*% x - b > 0 

# with solve.QP in version 2014.1-1 next does not work with EPS=1e-7
# [1] "Failure:  Error in projection"

ans <- spg(par=par0, fn=gaussmix.mloglik, gr=gaussmix.grad,
  project="projectLinear", projectArgs=list(A=Amat, b=b, meq=meq))

# # Does work!
# require("quadprog")
# 
# projectLinearOld <- function (par, A, b, meq){
#   n <- length(par)
#   if (meq > 0 | any(b - c(A %*% par) > 0)) {
#     ans <- solve.QP(Dmat = diag(1, n), dvec = rep(0, n),
#        Amat = t(A), bvec = b - c(A %*% par), meq = meq, factorized = TRUE)
#        par <- par + ans$solution
#     }
#   par
#   }
# 
# ans <- spg(par=par0, fn=gaussmix.mloglik, gr=gaussmix.grad,
#    project="projectLinearOld", projectArgs=list(A=Amat, b=b, meq=meq))

if (0 != ans$convergence) stop("test did not converge!")

fuzz <- 5e-5
if(fuzz < max(abs(ans$par -
   c( 0.2103359277577284137,  0.2191738028962620377, 0.2174358494266191433,
      0.3530544199193904609,  7.0060291485783237064, 11.2527073428970716407,
     -0.0166017473519236673,  2.9360474287487265954, 2.0609328632879644339)))){

#above is Mint 3.11.0-12-generic #19-Ubuntu SMP  x86_64
# Windows . using R version 3.1.1 (2014-07-10)
#     * using platform: x86_64-w64-mingw32 (64-bit)
# 'i386' 
#  [1]  0.2103348878325739246  0.2191732399236242523  0.2174403704585151642
#  [4]  0.3530515017852868809  7.0060519130170799684 11.2527120792098980218
#  [7] -0.0165690026560232316  2.9360806878664402753  2.0609368544772692644
#  arch 'x64' 
#  [1]  0.2103355761431467408  0.2191734750015318922  0.2174374428015042882
#  [4]  0.3530535060538171899  7.0060405011166100309 11.2527109132279594661
#  [7] -0.0165922248714034798  2.9360604372491612146  2.0609297716047816351

   print(ans$par, digits=18)
   cat("difference:\n")
   print(ans$par -
   c( 0.2103359277577284137,  0.2191738028962620377, 0.2174358494266191433,
      0.3530544199193904609,  7.0060291485783237064, 11.2527073428970716407,
     -0.0166017473519236673,  2.9360474287487265954, 2.0609328632879644339), 
    digits=18)
   stop("converged to different parameter values!")
   }

if(fuzz < max(abs(ans$value - 1388.64728677794915))){
   print(ans$value, digits=18)
   stop("converged to different function value!")
   }

ans

# $par
# [1]  0.21033488  0.21917325  0.21744042  0.35305145  7.00605176 11.25271191 # -0.01656842  2.93608090  2.06093738
# 
# $value
# [1] 1388.647
# 
# $gradient
# [1] 0.0001082927
# 
# $fn.reduction
# [1] 64.40326
# 
# $iter
# [1] 945
# 
# $feval
# [1] 1069
# 
# $convergence
# [1] 0
# 
# $message
# [1] "Successful convergence"


##########################
# Gaussian mixture density using upper and lower

Amat <- matrix(0, 1, 9)
Amat[1, 1:4] <- 1  # corresponds to equality

b <- 1
meq <- 1

ans2 <- spg(par=par0, fn=gaussmix.mloglik, gr=gaussmix.grad, 
   lower=c(rep(0,4), rep(-Inf, 4), 0), upper=c(rep(1,4), rep(Inf, 4), Inf),
   project="projectLinear", projectArgs=list(A=Amat, b=b, meq=meq))

if(fuzz < max(abs(ans$par - ans2$par))){
   print(ans2$par, digits=18)
   cat("difference:\n")
   print(ans$par - ans2$par, digits=18)
   stop("converged to different parameter values with lower and upper!")
   }

if(fuzz < max(abs(ans$value - ans2$value))){
   print(ans2$value, digits=18)
   stop("converged to different function value with lower and upper!")
   }

ans2

##########################
