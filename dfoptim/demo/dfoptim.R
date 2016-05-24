
#######################################################################################################
 rosbkext <- function(x){
# Extended Rosenbrock function
 n <- length(x)
 sum (100*(x[1:(n-1)]^2 - x[2:n])^2 + (x[1:(n-1)] - 1)^2)
 }

np <- 10
set.seed(123)

p.0 <- rnorm(np)
xm1 <- nmk(fn=rosbkext, par=p.0) # maximum `fevals' is not sufficient to find correct minimum
xm2 <- nmk(fn=rosbkext, par=p.0, control=list(maxfeval=5000)) # finds the correct minimum 
ans.optim <- optim(fn=rosbkext, par=p.0, method="Nelder-Mead", control=list(maxit=5000))   # terminates with inferior estimates
ans.hj <- hjk(fn=rosbkext, par=p.0)   # Hooke-Jeeves algorithm
xmb <- nmkb(fn=rosbkext, par=p.0, lower=-2, upper=2)
 
#######################################################################################################
### A non-smooth problem
nsf <- function(x) {
	f1 <- x[1]^2 + x[2]^2
	f2 <- x[1]^2 + x[2]^2 + 10 * (-4*x[1] - x[2] + 4)
	f3 <- x[1]^2 + x[2]^2 + 10 * (-x[1] - 2*x[2] + 6)
	max(f1, f2, f3)
}

p0 <- rnorm(3)
xm3 <- nmk(fn=nsf, par=p0)
xm3.hj <- hjk(fn=nsf, par=p0)

#######################################################################################################
### Another non-smooth problem
rosen <- function(x) {
# Rosen JB & Suzuki S (1965), Construction of non-linear programming test problems, Comm. ACM, 8, p. 113
	f1 <- x[1]^2 + x[2]^2 + 2*x[3]^2 + x[4]^2 - 5*x[1] - 5*x[2] - 21*x[3] + 7*x[4]
	f2 <- f1 + 10 * (sum(x^2) + x[1] - x[2] + x[3] - x[4] - 8)
	f3 <- f1 + 10 * (sum(x^2) + x[2]^2 + x[4]^2 - x[1] - x[4] - 10)
	f4 <- f1 + 10 * (sum(x^2) + x[1]^2 - x[4]^2 + 2*x[1] - x[2] - x[4] - 5)
	max(f1, f2, f3, f4)
}
# Global minimum value is -44 @ (0, 1, 2, -1)

p0 <- rnorm(4)
xm4 <- nmk(fn=rosen, par=p0)
xm4.hj <- hjk(fn=rosen, par=p0)
xm4b <- nmkb(fn=rosen, par=p0, lower=-2, upper=3)

#######################################################################################################
### Non-smooth problem #3
hald <- function(x) {
#Hald J & Madsen K (1981), Combined LP and quasi-Newton methods for minimax optimization, Mathematical Programming, 20, p.42-62.
	i <- 1:21
	t <- -1 + (i - 1)/10
	f <- (x[1] + x[2] * t) / ( 1 + x[3]*t + x[4]*t^2 + x[5]*t^3) - exp(t)
	max(abs(f))
	}
# Correct solution:  x* = (
# Minimum value = 0.002

p0 <- runif(5)
xm5 <- nmk(fn=hald, par=p0)
xm5.hj <- hjk(fn=hald, par=p0)
xm5b <- nmkb(fn=hald, par=p0, lower=c(0,0,0,0,-2), upper=4)
