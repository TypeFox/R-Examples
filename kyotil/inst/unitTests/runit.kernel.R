library("RUnit")
library("kyotil")


test.kernel <- function() {

RNGkind("Mersenne-Twister", "Inversion")
tolerance=1e-3
if(file.exists("C:/_checkReproducibility")) tolerance=1e-6


set.seed(1)
X = cbind(x1=rnorm(n=5), x2=rnorm(n=5))
dim(X)
X2 = cbind(x1=rnorm(n=3), x2=rnorm(n=3))
dim(X2)

K = getK(X,"linear")
dim(K)
checkEqualsNumeric(mean(K), 0.03497235, tol=tolerance)

K = getK(X,"linear",X2=X2)
dim(K)
K1 = getK(X2,"l",X2=X)
dim(K1)
checkTrue(all(K==t(K1)))


# RBF kernel
K = getK(X,"rbf",para=1,X2=X2)
checkEqualsNumeric(mean(K), 0.1994038, tol=tolerance)
K1 = getK(X2,"r",para=1,X2=X)
checkTrue(all(K==t(K1)))



# IBS kernel for ternary data 
X <- as.matrix(expand.grid(0:3,0:2))
checkException(getK(expand.grid(0:3,0:2),kernel = 'ibs'))
checkException(getK(as.matrix(expand.grid(0:3,0:2)),kernel = 'ibs'))

X <- as.matrix(expand.grid(0:2,0:2))
K = getK(X,kernel = 'ibs')
checkEqualsNumeric(mean(K), 0.5555556, tol=tolerance)

# add weight
set.seed(2)
w = runif(ncol(X))
K = getK(X,kernel = 'ibs',para = w) 
checkEqualsNumeric(mean(K^3),  0.3186303, tol=tolerance)
checkException(getK(X,kernel = 'ibs',para = -1))

# IBS kernel for binary data via option 'h' for 'hamming similarity measure'
X <- as.matrix(expand.grid(0:1,0:1))
K=getK(X,kernel = 'h')
checkEqualsNumeric(mean(K^3),.3125, tol=tolerance)
checkException(getK(X,kernel = 'h',para = -1))

# add weight
checkEqualsNumeric(mean(getK(X,kernel = 'h',para = w[]) ^3), 0.3762837, tol=tolerance)


n = 200
n2 = 100


k <- 6 # number of covariates
w <- runif(k) # weights for hamming and ibs kernels

# binary data {1,2} for hamming kernel
X.bin <- matrix(sample.int(2,n*k,replace = TRUE),n,k) - 1 
X2.bin = matrix(sample.int(2,n2*k,replace = TRUE),n2,k) - 1

checkEqualsNumeric(
    getK(X.bin,kernel = 'hamming',X2=X2.bin)
    ,
    getK(X.bin,kernel = 'hamming',X2=X2.bin,C = TRUE),
tolerance=tolerance)

checkEqualsNumeric(
    getK(X.bin,kernel = 'hamming',X2=X2.bin,para = w)
    ,
    getK(X.bin,kernel = 'hamming',X2=X2.bin,para = w,C = TRUE),
tolerance=tolerance)


# ternary data must be in {0,1,2} for ibs kernel
X.tern <- matrix(sample.int(3,n*k,replace = TRUE),n,k) - 1
X2.tern <- matrix(sample.int(3,n2*k,replace = TRUE),n2,k) - 1

checkEqualsNumeric(
    getK(X.tern,kernel = 'ibs',X2=X2.tern)
    ,
    getK(X.tern,kernel = 'ibs',X2=X2.tern,C = TRUE),
tolerance=tolerance)

checkEqualsNumeric(
    getK(X.tern,kernel = 'ibs',X2=X2.tern,para = w)
    ,
    getK(X.tern,kernel = 'ibs',X2=X2.tern,para = w,C = TRUE),
tolerance=tolerance)


X = cbind(x1=rnorm(n=n), x2=rnorm(n=n))
X2 = cbind(x1=rnorm(n=n2), x2=rnorm(n=n2))

checkEqualsNumeric(getK(X,kernel = 'linear',X2=X2),getK(X,kernel = 'linear',X2=X2,C = TRUE))
checkEqualsNumeric(getK(X,kernel = 'eucl',X2=X2),getK(X,kernel = 'eucl',X2=X2,C = TRUE))
checkEqualsNumeric(getK(X,kernel = 'poly',para = .5,X2=X2),getK(X,kernel = 'poly',para = .5,X2=X2,C = TRUE))
checkEqualsNumeric(getK(X,kernel = 'rbf',para = .5,X2=X2),getK(X,kernel = 'rbf',para = .5,X2=X2,C = TRUE))



}
