require("numDeriv")

##### Example 0
set.seed(123)
f <- function(x) {
  n <- length(x)
  f <- rep(NA, n)
  vec <- 1:(n-1)
  f[vec] <- x[vec]^2 + (-1)^vec * x[vec]*exp(x[vec+1])
  f[n] <- x[n]*exp(x[n])
  f
  }

x0 <- runif(5)
ans1 <- jacobian(func=f, x=x0,  method="complex")
print(ans1, digits=18)
#max.diff1:  3.571277e-11 
ans2 <- jacobian(func=f, x=x0)

err <- max(abs(ans1 - ans2))
cat("max.diff1: ", err, "\n")
if (1e-10 < err ) stop("Example 0 jacobian test failed.")


###### Example 1
broydt <- function(x, h=0.5) {
        n <- length(x)
        f <- numeric(n)
        f[1] <- ((3 - h*x[1]) * x[1]) - 2*x[2] + 1
        tnm1 <- 2:(n-1)
        f[tnm1] <- ((3 - h*x[tnm1])*x[tnm1]) - x[tnm1-1] - 2*x[tnm1+1] + 1
        f[n] <- ((3 - h*x[n]) * x[n]) - x[n-1] + 1
        sum(f*f)
    }


set.seed(123)
p0 <- runif(10)

ans1 <- grad(func=broydt, x=p0, method="complex")
#print(ans1, digits=18)
ans2 <- grad(func=broydt, x=p0)

err <- max(abs(ans1 - ans2))
cat("max.diff1: ", err, "\n")
#max.diff1:  4.977583e-10 
##max.diff1:  9.386859e-09 
if (1e-8 < err ) stop("broydt gradient test failed.")


h1 <- hessian(func=broydt, x=p0, method="complex")
#print(h1, digits=18)
h2 <- hessian(func=broydt, x=p0)
#print(h2, digits=18)
err <- max(abs(h1 - h2))
#print(err, digits=18)

cat("max.diff1: ", err , "\n")
#max.diff1:  9.386859e-09 
##max.diff1:  8.897979e-08 
if (1e-7 < err ) stop("broydt hessian test failed.")


###### Example 2
sc2.f <- function(x){
n <- length(x)
vec <- 1:n
sum(vec * (exp(x) - x)) / n
}

sc2.g <- function(x){
n <- length(x)
vec <- 1:n
vec * (exp(x) - 1) / n
}

sc2.h <- function(x){
n <- length(x)
hess <- matrix(0, n, n)
vec <- 1:n
diag(hess) <- vec*exp(x)/n
hess
}

set.seed(123)
#x0 <- rexp(10, rate=0.1)
x0 <- rnorm(100)
exact <- sc2.g(x0)

ans1 <- grad(func=sc2.f, x=x0, method="complex")
#print(ans1, digits=18)
err <- max(abs(exact - ans1)/(1 + abs(exact)))
err
#[1] 0
if (1e-14 < err ) stop("sc2 grad complex test failed.")


ans2 <- grad(func=sc2.f, x=x0)
err <- max(abs(exact - ans2)/(1 + abs(exact)))
err
# [1] 9.968372e-08
##[1] 9.968372e-08
if (1e-7 < err ) stop("sc2 grad Richardson test failed.")


exact <- sc2.h(x0)

system.time(ah1 <- hessian(func=sc2.f, x=x0, method="complex"))
#elapsed 4.14 
err <- max(abs(exact - ah1)/(1 + abs(exact)))
err
#  [1] 1.13183e-13
## [1] 1.13183e-13
if (1e-12 < err ) stop("sc2 hessian complex test failed.")

system.time(ah2 <- hessian(func=sc2.f, x=x0))
#elapsed  2.537 
err <- max(abs(exact - ah2)/(1 + abs(exact)))
err
# [1] 3.415308e-06
##[1] 6.969096e-08
if (1e-5 < err ) stop("sc2 hessian Richardson test failed.")


###### Example 3
rosbkext.f <- function(p, cons=10){
n <- length(p)
j <- 1: (n/2)
tjm1 <- 2*j - 1
tj <- 2*j 
sum (cons^2*(p[tjm1]^2 - p[tj])^2 + (p[tj] - 1)^2)
}

rosbkext.g <- function(p, cons=10){
n <- length(p)
g <- rep(NA, n)
j <- 1: (n/2)
tjm1 <- 2*j - 1
tj <- 2*j 
g[tjm1] <- 4*cons^2 * p[tjm1] * (p[tjm1]^2 - p[tj])
g[tj] <- -2*cons^2 * (p[tjm1]^2 - p[tj]) + 2 * (p[tj] - 1)
g
}

set.seed(123)
p0 <- runif(10)
exact <- rosbkext.g(p0, cons=10)

numd1 <- grad(func=rosbkext.f, x=p0, cons=10, method="complex") # not as good 
#print(numd1, digits=18)

err <- max(abs(exact - numd1)/(1 + abs(exact)))
err
# [1] 1.203382e-16
##[1] 1.691132e-16
if (1e-15 < err ) stop("rosbkext grad complex test failed.")


numd2 <- grad(func=rosbkext.f, x=p0, cons=10)
err <- max(abs(exact - numd2)/(1 + abs(exact)))
err
# [1] 5.825746e-11
##[1] 4.020598e-10
if (1e-9 < err ) stop("rosbkext grad Richardson test failed.")


###### Example 4
genrose.f <- function(x, gs=100){ 
  # objective function 
  ## One generalization of the Rosenbrock banana valley function (n parameters) 
  n <- length(x) 
  1.0 + sum (gs*(x[1:(n-1)]^2 - x[2:n])^2 + (x[2:n] - 1)^2) 
  }
 
genrose.g <- function(x, gs=100){ 
  # vectorized gradient for genrose.f # Ravi Varadhan 2009-04-03 
  n <- length(x)
  gg <- as.vector(rep(0, n)) 
  tn <- 2:n 
  tn1 <- tn - 1 
  z1 <- x[tn] - x[tn1]^2 
  z2 <- 1 - x[tn] 
  gg[tn] <- 2 * (gs * z1 - z2) 
  gg[tn1] <- gg[tn1] - 4 * gs * x[tn1] * z1 
  return(gg) 
  } 
  
#set.seed(123)
#p0 <- runif(10)
p0 <- rep(pi, 1000)
exact <- genrose.g(p0, gs=100)

numd1 <- grad(func=genrose.f, x=p0, gs=100, method="complex")
err <- max(abs(exact - numd1)/(1 + abs(exact)))
err
# [1] 2.556789e-16
##[1] 2.556789e-16
if (1e-15 < err ) stop("genrose grad complex test failed.")

numd2 <- grad(func=genrose.f, x=p0, gs=100) 
err <- max(abs(exact - numd2)/(1 + abs(exact)))
err
# [1] 1.847244e-09
##[1] 1.847244e-09
if (1e-8 < err ) stop("genrose grad Richardson test failed.")

##### Example 5
# function of single variable
fchirp <- function(x, b, k) exp(-b*x) * sin(k*x^4)
dchirp <- function(x, b, k) exp(-b*x) * (4 * k * x^3 * cos(k*x^4) - b * sin(k*x^4))

x <- seq(-3, 3, length=500)
y <- dchirp(x, b=1, k=4)
#plot(x, y, type="l")
y1 <- grad(func=fchirp, x=x, b=1, k=4, method="complex")
#lines(x, y1, col=2, lty=2)
err <- max(abs(y-y1))
err
# [1] 4.048388e-10
##[1] 4.048388e-10
if (1e-9 < err ) stop("chirp grad complex test failed.")

y2 <- grad(func=fchirp, x=x, b=1, k=4)
#lines(x, y2, col=3, lty=2)
err <- max(abs(y-y2))
err
# [1] 5.219681e-08
##[1] 5.219681e-08
if (1e-7 < err ) stop("chirp grad Richardson test failed.")

