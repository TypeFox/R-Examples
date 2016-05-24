## ----eval=FALSE----------------------------------------------------------
#  lm(y ~ x1 + x2 + ... + xm + f1 + f2 + ... + fn)

## ----eval=FALSE----------------------------------------------------------
#  felm(y ~ x1 + x2 + ... + xm | f1 + f2 + ... + fn)

## ----echo=FALSE----------------------------------------------------------
  cores <- as.integer(Sys.getenv('SG_RUN'))
  if(is.na(cores)) options(lfe.threads=1)

## ------------------------------------------------------------------------
library(lfe)
set.seed(42)
x1 <- rnorm(20)
f1 <- sample(8,length(x1),replace=TRUE)/10
f2 <- sample(8,length(x1),replace=TRUE)/10
e1 <- sin(f1) + 0.02*f2^2 + rnorm(length(x1))
y <-  2.5*x1 + (e1-mean(e1))
summary(est <- felm(y ~ x1 | f1 + f2))

## ------------------------------------------------------------------------
ef <- efactory(est)
is.estimable(ef,est$fe)
getfe(est)

## ------------------------------------------------------------------------
data.frame(f1,f2,comp=est$cfactor)

## ------------------------------------------------------------------------
f1 <- factor(f1); f2 <- factor(f2)
summary(lm(y ~ x1 + f1 + f2))

## ----eval=FALSE----------------------------------------------------------
#  est <- felm(logwage ~ x1 + x2 | id + firm + edu)
#  getfe(est)

## ----eval=FALSE----------------------------------------------------------
#  logwage ~ x1 + x2 + edu | id + firm

## ----eval=FALSE----------------------------------------------------------
#  logwage ~ x1 + x2 | firm + edu + id

## ----eval=FALSE----------------------------------------------------------
#  y ~ x1 + x2 + f1 + f2 + f3

## ----eval=FALSE----------------------------------------------------------
#  est <- felm(y ~ x1 + x2 | f1 + f2 + f3)

## ------------------------------------------------------------------------
library(lfe)
x1 <- rnorm(100)
f1 <- sample(7,100,replace=TRUE)
f2 <- sample(8,100,replace=TRUE)/8
f3 <- sample(10,100,replace=TRUE)/10
e1 <- sin(f1) + 0.02*f2^2  + 0.17*f3^3 + rnorm(100)
y <-  2.5*x1 + (e1-mean(e1))
summary(est <- felm(y ~ x1 | f1 + f2 + f3))

## ------------------------------------------------------------------------
ef <- function(gamma,addnames) {
  ref2 <- gamma[[8]]
  ref3 <- gamma[[16]]
  gamma[1:7] <- gamma[1:7]+ref2+ref3
  gamma[8:15] <- gamma[8:15]-ref2
  gamma[16:25] <- gamma[16:25]-ref3
  if(addnames) {
    names(gamma) <- c(paste('f1',1:7,sep='.'),
                          paste('f2',1:8,sep='.'),
                          paste('f3',1:10,sep='.'))
  }
  gamma
}
is.estimable(ef,fe=est$fe)
getfe(est,ef=ef)

## ------------------------------------------------------------------------
getfe(est)

## ------------------------------------------------------------------------
efactory(est)

## ------------------------------------------------------------------------
f1 <- factor(f1); f2 <- factor(f2); f3 <- factor(f3)
ef <- function(gamma,addnames) {
  ref1 <- gamma[[1]]
  ref2 <- gamma[[8]]
  ref3 <- gamma[[16]]
  # put the intercept in the first coordinate
  gamma[[1]] <- ref1+ref2+ref3
  gamma[2:7] <- gamma[2:7]-ref1
  gamma[8:14] <- gamma[9:15]-ref2
  gamma[15:23] <- gamma[17:25]-ref3
  length(gamma) <- 23
  if(addnames) {
    names(gamma) <- c('(Intercept)',paste('f1',levels(f1)[2:7],sep=''),
                          paste('f2',levels(f2)[2:8],sep=''),
                          paste('f3',levels(f3)[2:10],sep=''))
  }
  gamma
}
getfe(est,ef=ef,bN=1000,se=TRUE)
#compare with lm
summary(lm(y ~ x1 + f1 + f2 + f3))

## ------------------------------------------------------------------------
set.seed(128)
x1 <- rnorm(25)
f1 <- sample(9,length(x1),replace=TRUE)
f2 <- sample(8,length(x1),replace=TRUE)
f3 <- sample(8,length(x1),replace=TRUE)
e1 <- sin(f1) + 0.02*f2^2  + 0.17*f3^3 + rnorm(length(x1))
y <-  2.5*x1 + (e1-mean(e1))
summary(est <- felm(y ~ x1 | f1 + f2 + f3))

## ------------------------------------------------------------------------
ef <- efactory(est)
is.estimable(ef,est$fe)

## ------------------------------------------------------------------------
f1 <- factor(f1); f2 <- factor(f2); f3 <- factor(f3)
D <- makeDmatrix(list(f1,f2,f3))
dim(D)
ncol(D) - as.integer(rankMatrix(D))

## ------------------------------------------------------------------------
lfe:::rankDefic(list(f1,f2,f3))

## ------------------------------------------------------------------------
summary(est <- felm(y ~ x1 | f1 + f2 + f3, exactDOF=TRUE))

## ------------------------------------------------------------------------
summary(est <- felm(y ~ x1 +f3 | f1 + f2, exactDOF=TRUE))
getfe(est)

## ------------------------------------------------------------------------
summary(est <- felm(y ~ x1 | f1 + f3 + f2, exactDOF=TRUE))
is.estimable(efactory(est),est$fe)
getfe(est)


## ------------------------------------------------------------------------
data.frame(f1,f2,f3)[1,]

## ------------------------------------------------------------------------
summary(est <- lm(y ~ x1 + f1 + f2 + f3))

## ------------------------------------------------------------------------
fe <- list(f1,f2,f3)
wwcomp <- compfactor(fe, WW=TRUE)

## ------------------------------------------------------------------------
lfe:::rankDefic(fe)
nlevels(wwcomp)

## ------------------------------------------------------------------------
nlevels(interaction(compfactor(fe),wwcomp))
# pick the largest component:
wwdata <- data.frame(y, x1, f1, f2, f3)[wwcomp==1, ]
print(wwdata)

## ------------------------------------------------------------------------
set.seed(135)
x <- rnorm(10000)
f1 <- sample(1000,length(x),replace=TRUE)
f2 <- (f1 + sample(18,length(x), replace=TRUE)) %% 500
f3 <- (f2 + sample(9,length(x),replace=TRUE)) %% 500
y <- x + 1e-4*f1 + sin(f2^2) +
  cos(f3)^3 + 0.5*rnorm(length(x))
dataset <- data.frame(y,x,f1,f2,f3)
summary(est <- felm(y ~ x | f1 + f2 + f3,
             data=dataset, exactDOF=TRUE))

## ------------------------------------------------------------------------
nlevels(est$cfactor)
is.estimable(efactory(est), est$fe)
nrow(alpha <- getfe(est))

## ------------------------------------------------------------------------
lfe:::rankDefic(est$fe)

## ------------------------------------------------------------------------
wwcomp <- compfactor(est$fe,WW=TRUE)
nlevels(wwcomp)
wwset <- wwcomp == 1
sum(wwset)
summary(wwest <- felm(y ~ x | f1 + f2 + f3,
             data=dataset, subset=wwset, exactDOF=TRUE))

## ------------------------------------------------------------------------
lfe:::rankDefic(wwest$fe)

## ------------------------------------------------------------------------
nrow(wwalpha <- getfe(wwest))

## ------------------------------------------------------------------------
head(alpha)
head(wwalpha)

## ------------------------------------------------------------------------
table(wwalpha[,'fe'])

