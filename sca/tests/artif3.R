#### Artificial example with 3 or 4 true components :
library(sca)

Sig <- function(p, rho){ r <- diag(p); r[col(r) != row(r)] <- rho; r}
library(MASS)## -> mvrnorm() :
rmvN <- function(n,p, rho) mvrnorm(n, mu=rep(0,p), Sigma= Sig(p, rho))

set.seed(253)
## random matrix
mr <- matrix(rnorm(1000), 50, 20)

.proctime00 <- proc.time()

(scr0 <- sca(cor(mr)))
(scr <-  sca(cor(mr), q = 5, corblocks = 0.12))


##- Nr. 1 --- p = 3+2+4 = 9 ------------------

m3b <- cbind(rmvN(100, 3, 0.7),
             rmvN(100, 2, 0.9),
             rmvN(100, 4, 0.8))
## Show near block-structure of cor. matrix :
symnum(cor(m3b), lower.tri = FALSE)
sc3b <- sca(cor(m3b))
sc3b

sc3c.1 <- sca(cor(m3b), corblocks = 0.1)
## -> gives the 3 "true" block components
sc3c.1
str(sc3c.1)


##- Nr. 2 --- p = 12+6+2+10 = 30 ------------------

m4b <- cbind(rmvN(500, 12, 0.7),
             rmvN(500,  6, 0.9),
             rmvN(500,  2, 0.9),
             rmvN(500, 10, 0.8))
C4 <- cor(m4b)
## Show near block-structure of cor. matrix :
symnum(C4, lower.tri = FALSE)
sc4b <- sca(C4)
sc4b

sc4c.1 <- sca(C4, corblocks = 0.1)
## -> gives the 4 "true" block components
sc4c.1
str(sc4c.1)

cat('Time elapsed: ',proc.time() - .proctime00,'\n')
