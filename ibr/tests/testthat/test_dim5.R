context("testing basic modeling in dimension 5")
#### dataset
set.seed(46)
n <- 50
f<-function(x1,x2,x3,x4,x5){10*sin(pi*x1*x2)+20*(x3-0.5)^2+10*x4+5*x5}
x1<-runif(n)
x2<-runif(n)
x3<-runif(n)
x4<-runif(n)
x5<-runif(n)
Z <- f(x1,x2,x3,x4,x5)
sigma2 <- 0.1*var(Z)
y <- Z + rnorm(n,0,sqrt(sigma2))
don <-cbind.data.frame(y,x1,x2,x3,x4,x5)
x1<-runif(n)
x2<-runif(n)
x3<-runif(n)
x4<-runif(n)
x5<-runif(n)
dongrid <-cbind.data.frame(x1,x2,x3,x4,x5)
### expectations
## Kernel
res2 <- ibr(y~x1+x2+x3+x4+x5,data=don)
pres2 <- predict(res2,dongrid)
## loading data
load("resk.rda")
load("presk.rda")
test_that("modelling results with kernel default df", {
    expect_true(all(abs(res$fitted-res2$fitted)<1e-08))
    expect_true(all(abs(res$iter-res2$iter)<1e-08))
    expect_true(all(abs(res$initialdf-res2$initialdf)<1e-08))
    expect_true(all(abs(res$finaldf-res2$finaldf)<1e-08))
    expect_true(all(abs(res$bandwidth-res2$bandwidth)<1e-08))
    expect_true(all(abs(res$criteria-res2$criteria)<1e-03))
    expect_true(all(abs(pres-pres2)<1e-08))
})
## DS
res2 <- ibr(y~x1+x2+x3+x4+x5,data=don,smoother="ds")
pres2 <- predict(res2,dongrid)
## loading data
load("resds.rda")
load("presds.rda")
test_that("modelling results with DS splines default df", {
    expect_true(all(abs(res$fitted-res2$fitted)<1e-08))
    expect_true(all(abs(res$iter-res2$iter)<1e-08))
    expect_true(all(abs(res$initialdf-res2$initialdf)<1e-08))
    expect_true(all(abs(res$finaldf-res2$finaldf)<1e-08))
    expect_true(all(abs(res$bandwidth-res2$bandwidth)<1e-08))
    expect_true(all(abs(res$criteria-res2$criteria)<1e-03))
    expect_true(all(abs(pres-pres2)<1e-08))
})
## LOWRANK DS
res2 <- ibr(y~x1+x2+x3+x4+x5,data=don,smoother="lrds",rank=40)
## loading data
load("reslrds.rda")
test_that("modelling results with lowrank DS splines default df and rank=40", {
    expect_true(all(abs(res$fitted-res2$fitted)<1e-06))
    expect_true(all(abs(res$iter-res2$iter)<1e-06))
})
