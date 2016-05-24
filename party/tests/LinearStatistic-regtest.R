
set.seed(290875)
library("party")

### get rid of the NAMESPACE
attach(asNamespace("party"))

### 
###
###    Regression tests for linear statistics, expectations and covariances
###    
###    functions defined in file `./src/LinearStatistics.c'    

### tests for function C_LinearStatistic
### Linear Statistics
x = matrix(c(rep.int(1,4), rep.int(0,6)), ncol = 1)
y = matrix(1:10, ncol = 1)
weights = rep(1, 10)
linstat = LinearStatistic(x, y, weights)
stopifnot(isequal(linstat, sum(1:4)))

weights[1] = 0
linstat = LinearStatistic(x, y, weights) 
stopifnot(isequal(linstat, sum(2:4)))

xf <- gl(3, 10)
yf <- gl(3, 10)[sample(1:30)]
x <- sapply(levels(xf), function(l) as.numeric(xf == l))
colnames(x) <- NULL
y <- sapply(levels(yf), function(l) as.numeric(yf == l))
colnames(y) <- NULL
weights <- sample(1:30)
linstat <- LinearStatistic(x, y, weights) 
stopifnot(isequal(linstat, as.vector(t(x) %*% diag(weights) %*% y)))

xf <- factor(cut(rnorm(6000), breaks = c(-Inf, -2, 0.5, Inf)))
x <- sapply(levels(xf), function(l) as.numeric(xf == l))
yf <- factor(cut(rnorm(6000), breaks = c(-Inf, -0.5, 1.5, Inf)))
y <- sapply(levels(yf), function(l) as.numeric(yf == l))
weights <- rep(1, nrow(x))
colnames(x) <- NULL
colnames(y) <- NULL
weights <- rep(1, 6000)
linstat <- LinearStatistic(x, y, weights)
stopifnot(isequal(as.vector(table(xf, yf)), linstat))
stopifnot(isequal(as.vector(t(x)%*%y), linstat))


### tests for function C_ExpectCovarInfluence
eci <- ExpectCovarInfluence(y, weights)
isequal(eci@sumweights, sum(weights))
isequal(eci@expectation, drop(weights %*% y / sum(weights)))
ys <- t(t(y) - eci@expectation)
stopifnot(isequal(eci@covariance, (t(ys) %*% (weights * ys)) /
                  sum(weights)))

### tests for function C_ExpectCovarLinearStatistic
### Conditional Expectation and Variance (via Kruskal-Wallis statistic)

### case 1: p > 1, q = 1
group <- gl(3, 5)
x <- sapply(levels(group), function(l) as.numeric(group == l))
y <- matrix(1:15, ncol = 1)
weights <- rep(1, 15)

linstat <- LinearStatistic(x, y, weights)
expcov <- ExpectCovarLinearStatistic(x, y, weights)
KW <- quadformTestStatistic(linstat, expcov@expectation, expcov@covariance)
kts <- kruskal.test(y ~ group)$statistic
stopifnot(isequal(KW, kts))

### case 2: p = 1, q > 1
linstat <- LinearStatistic(y, x, weights)
expcov <- ExpectCovarLinearStatistic(y, x, weights)
KW <- quadformTestStatistic(linstat, expcov@expectation, expcov@covariance)
kts <- kruskal.test(y ~ group)$statistic
stopifnot(isequal(KW, kts))

### case 3: p = 1, q = 1
x <- x[,1,drop = FALSE]
linstat <- LinearStatistic(x, y, weights)
expcov <- ExpectCovarLinearStatistic(x, y, weights)
KW <- quadformTestStatistic(linstat, expcov@expectation, expcov@covariance)
kts <- kruskal.test(y ~ as.factor(x))$statistic
stopifnot(isequal(KW, kts))

### case 4: p > 1, q > 1 via chisq.test
n <- 900
xf <- gl(3, n / 3)
yf <- gl(3, n / 3)[sample(1:n)]
x <- sapply(levels(xf), function(l) as.numeric(xf == l))
colnames(x) <- NULL
y <- sapply(levels(yf), function(l) as.numeric(yf == l))
colnames(y) <- NULL
weights <- rep(1, n)
linstat <- LinearStatistic(x, y, weights)
expcov <- ExpectCovarLinearStatistic(x, y, weights)
chi <- quadformTestStatistic(linstat, expcov@expectation, expcov@covariance)
chis <- chisq.test(table(xf, yf))$statistic
stopifnot(isequal(round(chi, 1), round(chis, 1)))

### tests for function C_PermutedLinearStatistic
### Linear Statistics with permuted indices
x <- matrix(rnorm(100), ncol = 2)
y <- matrix(rnorm(100), ncol = 2)
weights <- rep(1, 50)
indx <- 1:50
perm <- 1:50
stopifnot(isequal(LinearStatistic(x, y, weights), 
                  PermutedLinearStatistic(x, y, indx, perm)))
x <- matrix(1:10000, ncol = 2)
y <- matrix(1:10000, ncol = 2)

for (i in 1:100) {
  indx <- sample(1:ncol(y), replace = TRUE)
  perm <- sample(1:ncol(y), replace = TRUE)

  stopifnot(isequal(as.vector(t(x[indx,]) %*% y[perm, ]),
                    PermutedLinearStatistic(x, y, indx, perm)))
}
