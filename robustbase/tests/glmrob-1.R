library(robustbase)

source(system.file("xtraR/ex-funs.R", package = "robustbase"))
source(system.file("test-tools-1.R",  package = "Matrix", mustWork=TRUE))# assert.EQ


###>> 1 ------------------- family = poisson ------------------------------------

### very simple model [with outliers]
set.seed(113)
y <- rpois(17, lambda = 4) ## -> target:  beta_0 = log(E[Y]) = log(4) = 1.386294

y[1:2] <- 99:100 # outliers
y
rm1 <- glmrob(y ~ 1, family = poisson, trace = TRUE,
              acc = 1e-13) # default is just 1e-4
## and check the robustness weights:
assert.EQ(c(0.0287933850640724, 0.0284930623638766,
		      0.950239140568007, 0.874115394740014),
		    local({w <- rm1$w.r; w[ w != 1 ] }), tol = 1e-14)
assert.EQ(coef(rm1), c("(Intercept)" = 1.41710946076738),tol = 1e-14)

cm1 <- glm   (y ~ 1, family = poisson, trace = TRUE)

rmMT <- glmrob(y ~ 1, family = poisson, trace = TRUE, method="MT")
(sMT <- summary(rmMT))

if(FALSE) # for manual digging:
debug(robustbase:::glmrobMqle)

allresid <- function(obj, types = c("deviance", "pearson", "working", "response"))
{
    sapply(types, residuals, object = obj)
}

okFit <- function(obj, check.attr=FALSE, ...) {
  all.equal(obj$y, obj$fitted.values + residuals(obj, "response"),
            check.attributes=check.attr, ...)
}

## check validity of several methods simultaneously:
y. <- model.response(model.frame(rm1))
stopifnot(okFit(cm1), okFit(rm1), y. == y)

alr.c <- allresid(cm1)
alr.r <- allresid(rm1)

## MM --- just for now --
plot(resid(cm1),                resid(rm1), asp=1); abline(0,1, col=2)
plot(resid(cm1,type="pearson"), resid(rm1, type="pearson"), asp=1); abline(0,1, col=2)
plot(resid(cm1,type="working"), resid(rm1, type="working"), asp=1); abline(0,1, col=2)

## leave away the outliers --
cm0 <- glm   (y ~ 1, family = poisson, trace = TRUE, subset = -(1:2))
plot(resid(cm0),                resid(rm1)[-(1:2)], asp=1); abline(0,1, col=2)
plot(resid(cm0,type="pearson"), resid(rm1, type="pearson")[-(1:2)], asp=1); abline(0,1, col=2)
plot(resid(cm0,type="working"), resid(rm1, type="working")[-(1:2)], asp=1); abline(0,1, col=2)
plot(resid(cm0,type="response"), resid(rm1, type="response")[-(1:2)], asp=1); abline(0,1, col=2)


## Use weights (slow convergence !)
w2 <- c(rep(1,8), rep(10,9))
rm2 <- glmrob(y ~ 1, family = poisson, trace = TRUE,
              weights = w2, maxit = 500, acc = 1e-10) # default is just 1e-4
## slow convergence
stopifnot(okFit(rm2))


###>> 2 ------------------- family = binomial -----------------------------------

## Using  *factor*  y ...
x <- seq(0,5, length = 120)
summary(px <- plogis(-5 + 2*x))
set.seed(7)
(f <- factor(rbinom(length(x), 1, prob=px)))

summary(m.c0 <- glm   (f ~ x, family = binomial))
summary(m.r0 <- glmrob(f ~ x, family = binomial))

## add outliers --- in y:
f. <- f
f.[i1 <- 2:3] <- 1
f.[i0 <- 110+c(1,7)] <- 0
        m.c1 <- glm   (f. ~ x, family = binomial)
summary(m.r1 <- glmrob(f. ~ x, family = binomial)) ## hmm, not so robust?
stopifnot(m.r1$w.r[c(i0,i1)] < 1/3, # well, at least down weighted
	  ## and coefficients change less :
	  (coef(m.r1) - coef(m.c0)) / (coef(m.c1) - coef(m.c0)) < 1)
assert.EQ(c("(Intercept)" = -3.10817337603974, x = 1.31618564057790),
	  coef(m.r1), tol= 1e-14, giveRE=TRUE)

y <- as.numeric(as.character(f.))
m.r2  <- BYlogreg(x0=x, y=y, trace=TRUE, maxhalf= 10)
m.r2A <- BYlogreg(x0=x, y=y, trace= 2  , maxhalf= 15)
## different.. but not so much:
iB <- 1:5
assert.EQ(m.r2A[iB], m.r2[iB], tol = .003, giveRE=TRUE)


assert.EQ(c("Intercept" = -2.9554950286, x = 1.2574679132),
          ## 32-bit{ada-5}  -2.95549502890363   1.25746791332613
	  m.r2$coef, tol=8e-10, giveRE=TRUE)# seen 5.316e-10   for --disable-long-double
assert.EQ( c(0.685919891749065, 0.256419206157062),
          ## 32-bit{ada-5}:
          ## 0.685919891858219, 0.256419206203016)
	  m.r2$sterror, tol=4e-9)# seen 1.025e-9   for --disable-long-double

data(foodstamp)
str(foodstamp)
## Model with 'income' instead of log(income+1)  is "interesting"
## because BYlogreg() needs  maxhalf > 10 for convergence!
m.fs0   <- glm   (participation ~ ., family=binomial, data=foodstamp)
m.fs0QL <- glmrob(participation ~ ., family=binomial, data=foodstamp)
y.fs <- foodstamp[,"participation"]
X.fs0 <- model.matrix(m.fs0)
head(X.fs0)
## (former default) maxhalf = 10  leads to too early convergence:
m.fsWBY. <- BYlogreg(x0=X.fs0, y=y.fs,
                     addIntercept=FALSE, trace=TRUE, maxhalf=10)
m.fs.BY. <- BYlogreg(x0=X.fs0, y=y.fs, initwml=FALSE,
                     addIntercept=FALSE, trace=TRUE, maxhalf=10)
m.fsWBY <- BYlogreg(x0=X.fs0, y=y.fs,
		    addIntercept=FALSE, trace=TRUE, maxhalf=18)
m.fs.BY <- BYlogreg(x0=X.fs0, y=y.fs, initwml=FALSE,
		    addIntercept=FALSE, trace=TRUE, maxhalf=18)

assert.EQ(m.fsWBY.[iB], m.fsWBY[iB], tol= 0.07)## almost 7% different
assert.EQ(m.fs.BY.[iB], m.fs.BY[iB], tol= 0.08)

foodSt <- within(foodstamp, { logInc <- log(1 + income) ; rm(income) })

m.fsML <- glm   (participation ~ ., family=binomial, data=foodSt)
m.fsQL <- glmrob(participation ~ ., family=binomial, data=foodSt)
X.fs <- model.matrix(m.fsML)
stopifnot(dim(X.fs) == c(150, 4)) # including intercept!
try(## FIXME -- Mahalanobis fails with singular matrix, here:
m.fsWBY <- BYlogreg(x0=X.fs, y=y.fs,
		    addIntercept=FALSE, trace=TRUE, maxhalf=18)
)
## maxhalf=18 is too much --> no convergence (in 1000 steps)
m.fs.BY <- BYlogreg(x0=X.fs, y=y.fs, initwml=FALSE,
		    addIntercept=FALSE, trace=TRUE, maxhalf=18)
signif(
    rbind(ML = coef(m.fsML),   QL =coef(m.fsQL),
          WBY0=coef(m.fsWBY.), BY0=coef(m.fs.BY.),
          WBY =coef(m.fsWBY ), BY =coef(m.fs.BY)
          )
    , 4)


if(FALSE) {
## *scaling* of X (  ?? <==> ??   'sigma1' ) ------------------

## no "W" (Mahalanobis fail because of *singular* X):
m.fs.BY100 <- BYlogreg(x0=100*X.fs, initwml=FALSE,
                       y=y.fs,
                       addIntercept=FALSE, trace=TRUE, maxhalf=18)
## ==> no convergence

X1c <- cbind(1, 100*X.fs[,-1])
m.fsWBY1c <- BYlogreg(x0=X1c, y=y.fs,
                      addIntercept=FALSE, trace=TRUE, maxhalf=18)
## ==> illegal singularity$kind

}## not yet

###-------- Gamma ------------

## Realistic "data" {from help(glmrob)}:
mu <- c(122.131, 53.0979, 39.9039, 33.9232, 28.007,
        24.923, 21.5747, 19.6971, 18.4516)
ns.resid <- c(-0.0338228, 0.0923228, 0.0525284, 0.0317426, -0.035954,
              0.00308925, -0.026637, -0.0353932, -0.0244761)
Vmu <- c(14915.9, 2819.38, 1592.32, 1150.78, 784.39,
         621.156, 465.467, 387.978, 340.462)
Hp2  <- robustbase:::Huberprop2
## Hp2. <- robustbase:::Huberprop2.

## was: phis <- 2^(-70:-1)  -- but that was *not* reliable (on 32-bit e.g.)
phis <- 2^(-42:-1)
H1 <- sapply(phis, function(phi)
    Hp2(phi, ns.resid=ns.resid, mu=mu, Vmu=Vmu, tcc = 1.345))
## H2 <- sapply(phis, function(phi)
##     Hp2.(phi, ns.resid=ns.resid, mu=mu, Vmu=Vmu, tcc = 1.345))
dput(signif(H1))
H2 <- c(9.91741,
        9.88674, 9.89438, 9.88674, 9.88961, 9.88961, 9.88961, 9.88984,
        9.88973, 9.88964, 9.8897, 9.88975, 9.88976, 9.88975, 9.88974,
        9.88974, 9.88974, 9.88974, 9.88974, 9.88974, 9.88974, 9.88974,
        9.88975, 9.88975, 9.88975, 9.33161, 8.70618, 8.39347, 8.23714,
        8.15902, 8.12006, 7.16275, 3.38703, -0.0879886, -2.3322, -4.16929,
        -5.26821, -5.80526, -6.04822, -6.11538, -6.02613, -5.66718)
          all.equal(H1,H2, tolerance = 0) # -> see 8.869e-7
stopifnot(all.equal(H1,H2, tolerance = 1e-5))

if(dev.interactive(TRUE)) # shows that phi < 1e-12 is doubtful
  matplot(phis, cbind(H1,H2), log="x", ylim = rrange(H1), type="o")

