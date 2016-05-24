library(robustbase)

### Binomial example with  *small*  ni

N <- 51
set.seed(123)
table(ni <- rpois(N, lam=4))# has 4 '1's, (no '0')
n0 <- ni; n0[print(which(ni == 1)[1:2])] <- 0 # has two '0's
x <- seq(0,1, length=N)
pr.x <- plogis(5*(x - 1/2))
k  <- rbinom(N, size = ni, prob = pr.x)
k0 <- rbinom(N, size = n0, prob = pr.x)
cbind(k,ni, k0,n0)
g1 <- glm(cbind(k , ni-k ) ~ x, family = binomial)
coef(summary(g1))[,1:2]
g0  <- glm(cbind(k0, n0-k0) ~ x, family = binomial)# works too
g0. <- glm(cbind(k0, n0-k0) ~ x, family = binomial, subset = n0 > 0)
## all.equal(g0, g0.)
stopifnot(all.equal(print(coef(summary(g0))), coef(summary(g0.))))


rg1  <- glmrob(cbind(k , ni-k ) ~ x, family = binomial)
rg1. <- glmrob(cbind(k , ni-k ) ~ x, family = binomial,
               acc = 1e-10) # default is just 1e-4

stopifnot(all.equal(unname(coef(rg1.)), c(-2.37585864, 4.902389143), tolerance=1e-9),
          all.equal(coef(rg1),  coef(rg1.), tolerance=1e-4),
          all.equal(vcov(rg1.), vcov(rg1), tolerance = 1e-4))
rg1$iter
which(rg1.$w.r != 1) ## 7 of them :
str(rg1.["family" != names(rg1.)])

rg2 <- glmrob(cbind(k , ni-k ) ~ x, family = binomial,
               acc = 1e-10, tcc = 3) # large cutoff: almost classical
vcov(rg2) # << already close to limit
rg10 <- glmrob(cbind(k , ni-k ) ~ x, family = binomial, tcc = 10)
rgL  <- glmrob(cbind(k , ni-k ) ~ x, family = binomial, tcc = 100)

no.comp <- - match(c("call", "data", "family", "control", "tcc"), names(rg10))
stopifnot(all.equal(rg10[no.comp], rgL[no.comp], tolerance= 1e-14))

vcov(rgL) # is now the same as the following:
if(FALSE) { ## tcc=Inf fails: non-convergence / singular matrix from GOTO/Atlas3
 rgI <- glmrob(cbind(k , ni-k ) ~ x, family = binomial, tcc = Inf)
 ## tcc = Inf  still *FAILS* (!)
 stopifnot(all.equal(rgL[no.comp], rgI[no.comp], tolerance= 0))
 ## and is quite close to the classic one:
 (all.equal(vcov(rgI), vcov(g1)))
}

rg0 <-  glmrob(cbind(k0, n0-k0) ~ x, family = binomial)
## --> warning..
rg0. <- glmrob(cbind(k0, n0-k0) ~ x, family = binomial, subset = n0 > 0)

coef(summary(rg0)) # not yet good (cf. 'g0' above!) -- but the one of rg0. is
stopifnot(all.equal(coef(rg0), coef(rg0.)))


### Example where all ni >= 3  -- works better, now also correct as.var. !!
### ----------------- =======

min(n3 <- ni + 2)# = 3
k3 <- rbinom(N, size = n3, prob = pr.x)
g3 <- glm(cbind(k3 , n3-k3) ~ x, family = binomial)
(cfg <- coef(summary(g3))[,1:2])
stopifnot(all.equal(sqrt(diag(vcov(g3))), cfg[,2]))

rg3 <- glmrob(cbind(k3 , n3-k3) ~ x, family = binomial)
(s3 <- summary(rg3))
summary(rg3$w.r)
rg3.5 <- glmrob(cbind(k3 , n3-k3) ~ x, family = binomial, tcc = 5)
(s3.5 <- summary(rg3.5))
summary(rg3.5$w.r)# all 1
stopifnot(all.equal(coef(s3)[,1:2], coef(s3.5)[,1:2], tolerance = 0.02))

rg3.15 <- glmrob(cbind(k3 , n3-k3) ~ x, family = binomial, tcc = 15, acc=1e-10)
(s3.15 <- summary(rg3.15))

stopifnot(all.equal(coef(s3.15)[,1:2], cfg, tolerance = 1e-5),# 2e-6
          all.equal(cfg[,"Estimate"], rg3.15$coeff, tolerance= 1e-8) # 6.05e-10
          )
##rg3.15$eff # == 1

## doesn't change any more:
rg3.1000 <- glmrob(cbind(k3 , n3-k3) ~ x, family = binomial, tcc = 1000,
                   acc=1e-10)
stopifnot(all.equal(rg3.1000[no.comp],
                    rg3.15  [no.comp], tol = 1e-13))

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
