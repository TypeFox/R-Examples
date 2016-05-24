context("impacts")

test_that("optimized total impacts calculation is correct", {  
# Setup 
require(spatialprobit)
n <- 10
rho <- 0.7
X <- cbind(1, seq(-1, 1, length.out=10))
beta <- c(-1, 1)
I_n <- sparseMatrix(i=1:n, j=1:n, x=1)
W <- bandSparse(10, k = c(-1,1), diag = list(c(rep(0.5, 8), 1), c(1, rep(0.5, 8))))
S <- (I_n - rho * W)
QR <- qr(S)
ones <- rep(1, n)

# theoretical total impacts, see LeSage (2009)
mu <- as.double(solve(QR, X %*% beta))            # n x 1
D  <- diag(dnorm(mu))                             # n x n 
# for all non-constant variables (r=1)
S_r <- D %*% solve(S) * beta[2]                   # S_r(W)

# total impacts = 1/n * 1'_n %*% S_r(W) %*% 1_n = 1/n * sum(S_r(W))
total_impacts1 <- 1/n * sum(S_r)                  # 0.2287129

# optimized way to calculate total impacts via QR decomposition of S
# allows for vectorized computation without having a full n x n matrix
total_impacts2 <- colMeans(D %*% solve(QR, ones)) * beta[2]

expect_equal(total_impacts1,total_impacts2)
}) 