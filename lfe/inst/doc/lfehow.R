## ----eval=FALSE----------------------------------------------------------
#  Y ~ X1 + X2 + ... + Xk + D1 + D2 + ... + De

## ----eval=FALSE----------------------------------------------------------
#  felm(Y ~ X1 + X2 + ... + Xk | D1 + D2 + ... + De)

## ----eval=FALSE----------------------------------------------------------
#  Y ~ X1 + X2 | X3:D1 + D2 + D3

## ------------------------------------------------------------------------
set.seed(41)
x <- rnorm(500)
x2 <- rnorm(length(x))
x3 <- rnorm(length(x))

## ------------------------------------------------------------------------
f1 <- factor(sample(7,length(x),replace=TRUE))
f2 <- factor(sample(4,length(x),replace=TRUE))
f3 <- factor(sample(3,length(x),replace=TRUE))
eff1 <- rnorm(nlevels(f1))
eff2 <- rexp(nlevels(f2))
eff3 <- runif(nlevels(f3))

## ------------------------------------------------------------------------
 y <- x + 0.5*x2 + 0.25*x3 + eff1[f1] + eff2[f2] + eff3[f3] + rnorm(length(x))

## ------------------------------------------------------------------------
demean <- function(v,fl) {
  Pv <- v; oldv <- v-1
  while(sqrt(sum((Pv-oldv)**2)) >= 1e-7) {
    oldv <- Pv
    for(f in fl) Pv <- Pv - ave(Pv,f)
  }
 Pv
}

## ------------------------------------------------------------------------
 fl <- list(f1,f2,f3)
 Py <- demean(y,fl)
 Px <- demean(x,fl)
 Px2 <- demean(x2,fl)
 Px3 <- demean(x3,fl)

## ------------------------------------------------------------------------
summary(lm(Py ~ Px + Px2 + Px3 - 1))

## ----echo=FALSE----------------------------------------------------------
  cores <- as.integer(Sys.getenv('SG_RUN'))
  if(is.na(cores)) options(lfe.threads=1)

## ------------------------------------------------------------------------
 library(lfe, quietly=TRUE)
 summary(est <- felm(y ~ x + x2 + x3 | f1+f2+f3))

## ----tidy=FALSE----------------------------------------------------------
ef <- function(v,addnames) {
  r1 <- v[[1]]
  r2 <- v[[8]]
  r3 <- v[[12]]
  result <- c(r1+r2+r3,v[2:7]-r1,v[9:11]-r2,v[13:14]-r3)
  if(addnames) names(result) <- c('(Intercept)',
                             paste('f1',2:7,sep='.'),
                             paste('f2',2:4,sep='.'),
                             paste('f3',2:3,sep='.'))
  result
}
# verify that it's estimable
is.estimable(ef,est$fe)
getfe(est, ef=ef, se=TRUE, bN=10000)

## ------------------------------------------------------------------------
summary(lm(y ~ x + x2 + x3 + f1 + f2 + f3))

