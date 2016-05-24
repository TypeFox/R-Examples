### test subsample
### LU decomposition and singular subsamples handling
require(robustbase)
source(system.file("xtraR/subsample-fns.R", package = "robustbase", mustWork=TRUE))
source(system.file("test-tools-1.R", package="Matrix", mustWork=TRUE))
require(Matrix)

cat("doExtras:", doExtras <- robustbase:::doExtras(),"\n")
showProc.time()

A <- matrix(c(0.001, 1, 1, 2), 2)
set.seed(11)
str(sa <- tstSubsample(A))

A <- matrix(c(3, 2, 6, 17, 4, 18, 10, -2, 12), 3)
tstSubsample(A)

## test some random matrix
set.seed(1002)
A <- matrix(rnorm(100), 10)
tstSubsample(A)

## test singular matrix handling
A <- matrix(c(1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1), 4, byrow=TRUE)
tstSubsample(A)


## test subsample with mts > 0
data <- data.frame(y = rnorm(9), expand.grid(A = letters[1:3], B = letters[1:3]))
x <- model.matrix(y ~ ., data)
y <- data$y
## this should produce a warning and return status == 2
showSys.time(z <- Rsubsample(x, y, mts=2))
stopifnot(z$status == 2)


## test equilibration
## columns only
X <- matrix(c(1e-7, 2, 1e-10, 0.2), 2)
y <- 1:2
tstSubsample(t(X), y)

## rows only
X <- matrix(c(1e-7, 2, 1e10, 0.2), 2)
y <- 1:2
tstSubsample(X, y)

## both
X <- matrix(c(1e-7, 1e10, 2, 2e12), 2)
y <- 1:2
tstSubsample(X, y)
showProc.time()


## test real data example
data(possumDiv)## 151 * 9; the last two variables are factors
with(possumDiv, table(eucalyptus, aspect))

mf <- model.frame(Diversity ~ .^2, possumDiv)
X <- model.matrix(mf, possumDiv)
y <- model.response(mf)
stopifnot(qr(X)$rank == ncol(X))

## this used to fail: different pivots in step 37
str(s1 <- tstSubsample(X, y))
s2 <- tstSubsample(X / max(abs(X)), y / max(abs(X)))
s3 <- tstSubsample(X * 2^-50, y * 2^-50)
## all components *BUT*  x, y, lu, Dr, Dc, rowequ, colequ :
nm <- names(s1); nm <- nm[is.na(match(nm, c("x","y","lu", "Dr", "Dc", "rowequ", "colequ")))]
stopifnot(all.equal(s1[nm], s2[nm], tolerance=1e-10),
	  all.equal(s1[nm], s3[nm], tolerance=1e-10))
showProc.time()

set.seed(10)
nsing <- sum(replicate(if(doExtras) 200 else 20, tstSubsampleSing(X, y)))
stopifnot(nsing == 0)
showProc.time()

## test example with many categorical predictors
set.seed(10)
r1 <- lmrob(Diversity ~ .^2 , data = possumDiv, cov="none")
## lmrob.S used to fail for this seed:
set.seed(108)
r2 <- lmrob(Diversity ~ .^2 , data = possumDiv, cov="none") #, trace=4)
showProc.time()

## investigate problematic subsample:
idc <- 1 + c(140, 60, 12, 13, 89, 90, 118, 80, 17, 134, 59, 94, 36,
         43, 46, 93, 107, 62, 57, 116, 11, 45, 35, 38, 120, 34, 29,
         33, 147, 105, 115, 92, 61, 91, 104, 141, 138, 129, 130, 84,
         119, 132, 6, 135, 112, 16, 67, 41, 102, 76, 111, 82, 148, 24,
         131, 10, 96, 0, 87, 21, 127, 56, 124)

rc <- lm(Diversity ~ .^2 , data = possumDiv, subset = idc)

X <- model.matrix(rc)
y <- possumDiv$Diversity[idc]
tstSubsample(X, y)## have different pivots ... could not find non-singular

lu <- LU.gaxpy(t(X))
stopifnot(lu$sing)
zc <- Rsubsample(X, y)
stopifnot(zc$status > 0)
## column 52 is linearly dependent and should have been discarded
## qr(t(X))$pivot

image(as(round(zc$lu -      (lu$L + lu$U - diag(nrow(lu$U))), 10), "Matrix"))
image(as( sign(zc$lu) - sign(lu$L + lu$U - diag(nrow(lu$U))),      "Matrix"))
showProc.time()

## test equilibration
## colequ only
X <- matrix(c(1e-7, 2, 1e-10, 0.2), 2)
y <- 1:2
tstSubsample(t(X), y)

## rowequ only
X <- matrix(c(1e-7, 2, 1e10, 0.2), 2)
y <- 1:2
tstSubsample(X, y)

## both
X <- matrix(c(1e-7, 1e10, 2, 2e12), 2)
y <- 1:2
tstSubsample(X, y)
showProc.time()

### real data, see MM's ~/R/MM/Pkg-ex/robustbase/hedlmrob.R
##  close to singular cov():
attach(system.file("external", "d1k27.rda", package="robustbase", mustWork=TRUE))

fm1 <- lmrob(y ~ a + I(a^2) + tf + I(tf^2) + A + I(A^2) + . , data = d1k27)
##     ^^^^^ gave error, earlier, now with a warning -- use ".vcov.w"
## --> cov = ".vcov.w"
fm2 <- lmrob(y ~ a + I(a^2) + tf + I(tf^2) + A + I(A^2) + . , data = d1k27,
             cov = ".vcov.w", trace = TRUE)
showProc.time()# 2.77

if(doExtras) {##-----------------------------------------------------------------

## Q: does it change to use numeric instead of binary factors ?
## A: not really ..
d1k.n <- d1k27
d1k.n[-(1:5)] <- lapply(d1k27[,-(1:5)], as.numeric)

fm1.n <- lmrob(y ~ a + I(a^2) + tf + I(tf^2) + A + I(A^2) + . , data = d1k.n)
fm2.n <- lmrob(y ~ a + I(a^2) + tf + I(tf^2) + A + I(A^2) + . , data = d1k.n,
             cov = ".vcov.w", trace = 2)

print(summary(weights(fm1, type="robustness")))
hist(weights(fm1, type="robustness"), main="robustness weights of fm1")
rug(weights(fm1, type="robustness"))
showProc.time()## 2.88

##
fmc <- lm   (y ~ poly(a,2)-a + poly(tf, 2)-tf + poly(A, 2)-A + . , data = d1k27)
print(summary(fmc))
## -> has NA's for  'a, tf, A'  --- bad that it did *not* work to remove them

nform <- update(formula(fm1), ~ .
                +poly(A,2)  -A  -I(A^2)
                +poly(a,2)  -a  -I(a^2)
                +poly(tf,2) -tf -I(tf^2))

fm1. <- lmrob(nform, data = d1k27)# now w/o warning !? !!
fm2. <- lmrob(nform, data = d1k27, cov = ".vcov.w", trace = TRUE)

## now lmrob takes care of NA coefficients automatically
print(lmrob(y ~ poly(a,2)-a + poly(tf, 2)-tf + poly(A, 2)-A + . , data = d1k27))
showProc.time() ## 4.24
} ## only if(doExtras) ##--------------------------------------------------------

## test exact fit property
set.seed(20)
data <- data.frame(y=c(rep.int(0, 20), rnorm(5)),
                   group=rep(letters[1:5], each=5))
x <- model.matrix(y ~ group, data)
lmrob.S(x, data$y, lmrob.control())
(ret <- lmrob(y ~ group, data))
summary(ret)


showProc.time()

