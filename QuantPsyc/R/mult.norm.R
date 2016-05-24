"mult.norm" <-
function (x,s=var(x), chicrit=.005)  
{ 
 
x <- as.matrix(x)
x <- na.omit(x)
n <- nrow(x)
p <- ncol(x)
dfchi <- p*(p+1)*(p+2)/6

# create 1 matrix
one <- matrix(1,n)

# create idendity matrix n*n
id <- matrix(0,n,n)
diag(id) <- 1

# create Q matrix ( I - 1/n * 1[n] * 1[n]' )
Q <- id - 1/n * one %*% t(one)

# create g matrix 
g <- Q %*% x %*% solve(s) %*% t(x) %*% Q


# mahalanobis distance = 
mahalanobis <- diag(g)
# Extreme D^2 =
ext <- qchisq(1-chicrit,p)

# b1 & b2 used for skew & kurt measure
b1p <- 1/(n^2) * sum(g^3)
b2p <- sum(diag(g)^2)/n

# Kappa measures
k1 <- n * b1p / 6
k2 <- (b2p - p*(p+2))/(8*p*(p+2)/n)^.5

# p-values
pskew <- 1 - pchisq(k1,dfchi)
pkurt <- 2 * ( 1 - pnorm(abs(k2)))

# create matrix to return
mat <- matrix(,2,3)
colnames(mat) <- c("Beta-hat", "kappa", "p-val")
rownames(mat) <- c("Skewness","Kurtosis")
mat[1,] <- c(as.numeric(b1p),as.numeric(k1),as.numeric(pskew))
mat[2,] <- c(as.numeric(b2p),as.numeric(k2),as.numeric(pkurt))
mult <- list(mult.test=mat,Dsq=mahalanobis, CriticalDsq=as.numeric(ext))
return(mult)

}

