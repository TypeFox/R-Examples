## Test implementation of M-S estimator
require(robustbase)
source(system.file("xtraR/m-s_fns.R", package = "robustbase", mustWork=TRUE))
source(system.file("xtraR/ex-funs.R", package = "robustbase", mustWork=TRUE))
source(system.file("test-tools-1.R",  package = "Matrix",     mustWork=TRUE))# assert.EQ

## dataset with factors and continuous variables:
data(education)
education <- within(education, Region <- factor(Region))
## for testing purposes:
education2 <- within(education, Group <- factor(rep(1:3, length.out=length(Region))))

## Test splitFrame (type fii is the only problematic type)
testFun <- function(formula, x1.idx) {
    obj <- lm(formula, education2)
    mf <- obj$model
    ret <- splitFrame(mf, type="fii")
    if (missing(x1.idx)) {
        print(ret$x1.idx)
        return(which(unname(ret$x1.idx)))
    }
    stopifnot(identical(x1.idx, which(unname(ret$x1.idx))))
}
testFun(Y ~ 1, integer(0))
testFun(Y ~ X1*X2*X3, integer(0))
testFun(Y ~ Region + X1 + X2 + X3, 1:4)
testFun(Y ~ 0 + Region + X1 + X2 + X3, 1:4)
testFun(Y ~ Region*X1 + X2 + X3, c(1:5, 8:10))
testFun(Y ~ Region*X1 + X2 + X3 + Region*Group, c(1:5, 8:18))
testFun(Y ~ Region*X1 + X2 + X3 + Region*Group*X2, c(1:6, 8:29))
testFun(Y ~ Region*X1 + X2 + Region*Group*X2, 1:28)
testFun(Y ~ Region*X1 + X2 + Region:Group:X2, 1:21)
testFun(Y ~ Region*X1 + X2*X3 + Region:Group:X2, c(1:6, 8:10, 12:23))
testFun(Y ~ (X1+X2+X3+Region)^2, c(1:7,10:12,14:19))
testFun(Y ~ (X1+X2+X3+Region)^3, c(1:19, 21:29))
testFun(Y ~ (X1+X2+X3+Region)^4, 1:32)
testFun(Y ~ Region:X1:X2 + X1*X2, c(1:1, 4:7))


control <- lmrob.control()
f.lm <- lm(Y ~ Region + X1 + X2 + X3, education)
splt <- splitFrame(f.lm$model)
y <- education$Y

## test orthogonalizing
x1 <- splt$x1
x2 <- splt$x2
tmp <- lmrob.lar(x1, y, control)
y.tilde <- tmp$resid
t1 <- tmp$coef
x2.tilde <- x2
T2 <- matrix(0, nrow=ncol(x1), ncol=ncol(x2))
for (i in 1:ncol(x2)) {
    tmp <- lmrob.lar(x1, x2[,i], control)
    x2.tilde[,i] <- tmp$resid
    T2[,i] <- tmp$coef
}

set.seed(10)
mss1 <- m_s_subsample(x1, x2.tilde, y.tilde, control, orth = FALSE)
mss1 <- within(mss1, b1 <- drop(t1 + b1 - T2 %*% b2))
set.seed(10)
mss2 <- m_s_subsample(x1, x2,       y,       control, orth = TRUE)
stopifnot(all.equal(mss1, mss2))

res <- vector("list", 100)
set.seed(0)
time <- system.time(for (i in seq_along(res)) {
    tmp <- m_s_subsample(x1, x2.tilde, y.tilde, control, FALSE)
    res[[i]] <- unlist(within(tmp, b1 <- drop(t1 + b1 - T2 %*% b2)))
})
cat('Time elapsed in subsampling: ', time,'\n')
## show a summary of the results  {"FIXME": output is platform dependent}
summary(res1 <- do.call(rbind, res))
## compare with fast S solution
fmS <- lmrob(Y ~ Region + X1 + X2 + X3, education, init="S")
coef(fmS)
fmS$scale

###  Comparing m-s_descent implementations()  {our C and R} : -------------------

ctrl <- control
#ctrl$trace.lev <- 5
ctrl$k.max <- 1
mC <- m_s_descent      (x1, x2, y, ctrl, mss2$b1, mss2$b2, mss2$scale+10)
mR <- m_s_descent_Ronly(x1, x2, y, ctrl, mss2$b1, mss2$b2, mss2$scale+10)
nm <- c("b1","b2", "scale", "res")
stopifnot(all.equal(mC[nm], mR[nm], check.attributes = FALSE, tolerance=5e-15))

## control$k.m_s <- 100
res3 <- vector("list", 100)
time <- system.time(for (i in seq_along(res3)) {
    ri <- res[[i]]
    res3[[i]] <- unlist(m_s_descent(x1, x2, y, control,
				    ri[1:4], ri[5:7], ri[8]))
})
cat('Time elapsed in descent proc: ', time,'\n')

## show a summary of the results   {"FIXME": output is platform dependent}
res4 <- do.call(rbind, res3)
summary(res4[,1:8])

stopifnot(all.equal( # 'test', not only plot:
	  res1[, "scale"],   res4[,"scale"], tol = 0.03),
	  res1[, "scale"] >= res4[,"scale"] - 1e-7 ) # 1e-7 just in case
     plot(res1[, "scale"],   res4[,"scale"])
abline(0,1, col=adjustcolor("gray", 0.5))

## Test lmrob.M.S
x <- model.matrix(fmS)
control$trace.lev <- 3
##      ---------   --
set.seed(1003)
fMS <- lmrob.M.S(x, y, control, fmS$model)
resid <- drop(y - x %*% fMS$coef)
assert.EQ(resid, fMS$resid, check.attributes=FALSE, tol = 1e-12)

## Test direct call to lmrob
## 1. trace_lev output:
set.seed(17)
fMS <- lmrob(Y ~ Region + X1 + X2 + X3, education, init = "M-S", trace.lev=2)

set.seed(13)
fiMS <- lmrob(Y ~ Region + X1 + X2 + X3, education, init = "M-S")
out2 <- capture.output(summary(fiMS))
writeLines(out2)

set.seed(13)
fiM.S <- lmrob(Y ~ Region + X1 + X2 + X3, education, init=lmrob.M.S)
out3 <- capture.output(summary(fiM.S))

## must be the same {apart from the "init=" in the call}:
i <- 3
stopifnot(identical(out2[-i], out3[-i]))
## the difference:
c(rbind(out2[i], out3[i]))


###  "Skipping design matrix equilibration" warning can arise for reasonable designs -----
set.seed(1)
x2 <- matrix(rnorm(2*30), 30, 2)
data <- data.frame(y = rnorm(30), group = rep(letters[1:3], each=10), x2)

obj <- lmrob(y ~ ., data, init="M-S", trace.lev=1)

## illustration: the zero row is introduced during the orthogonalization of x2 wrt x1
## l1 regression always produces p zero residuals
## by chance, the zero residuals of multiple columns happen to be on the same row
sf <- splitFrame(obj$model)
x1 <- sf$x1
x2 <- sf$x2
control <- obj$control

## orthogonalize
x2.tilde <- x2

for(i in 1:ncol(x2)) {
    tmp <- lmrob.lar(x1, x2[,i], control)
    x2.tilde[,i] <- tmp$resid
}
x2.tilde == 0


## Specifying init="M-S" for a model without categorical variables
## used to cause a segfault; now uses "S"
lmrob(LNOx ~ LNOxEm, NOxEmissions[1:10,], init="M-S")

## Now an ANOVA model with *only* categorical variables
n <- 64 # multiple of 16
stopifnot(n %% 16 == 0)
d.AOV <- data.frame(y = round(100*rnorm(64)),
		    A=gl(4,n/4), B=gl(2,8, n), C=gl(2,4,n))
fm <- lmrob(y ~ A*B*C, data = d.AOV, init = "M-S", trace.lev=2)

## lmrob_M_S(n = 64, nRes = 500, (p1,p2)=(16,0), (orth,subs,desc)=(1,1,1))
##  Starting subsampling procedure.. Error in lmrob.M.S(x, y, control, mf) :
##   'Calloc' could not allocate memory (18446744073709551616 of 4 bytes)

## BTW: Can we compute an  M-estimate (instead of MM-*) as we
## ---  cannot have any x-outliers in such an ANOVA!
