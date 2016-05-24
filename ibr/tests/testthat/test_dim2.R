context("testing basic modeling in dimension 2")
#### dataset
f <- function(x, y) { .75*exp(-((9*x-2)^2 + (9*y-2)^2)/4) +
                      .75*exp(-((9*x+1)^2/49 + (9*y+1)^2/10)) +
                      .50*exp(-((9*x-7)^2 + (9*y-3)^2)/4) -
                      .20*exp(-((9*x-4)^2 + (9*y-7)^2)) }
# define a (fine) x-y grid and calculate the function values on the grid
ngrid <- 5; xf <- seq(0,1, length=ngrid+2)[-c(1,ngrid+2)]
yf <- xf ; zf <- outer(xf, yf, f)
grid <- cbind(rep(xf, ngrid), rep(xf, rep(ngrid, ngrid)))
#generate a data set with function f and noise to signal ratio 5
noise <- .2 ; N <- 110 
set.seed(2557)
xr <- sort(runif(10,0.05,0.95)) ;
yr <- sort(runif(11,0.05,0.95))
zr <- outer(xr,yr,f)
set.seed(25)
std <- sqrt(noise*var(as.vector(zr))) ; noise <- rnorm(length(zr),0,std)
Z <- zr + matrix(noise,nrow(zr),ncol(zr))
# transpose the data to a column format 
xc <- rep(xr, length(yr))
yc <- rep(yr, each=length(xr))
X <- cbind(xc, yc) ; Zc <- as.vector(Z)
## data-frame
don <- cbind.data.frame(x=xc,y=yc,z=Zc)
## grid as data-frame
dongrid <- cbind.data.frame(x=grid[,1],y=grid[,2])

### expectations
## Kernel
res2 <- ibr(z~x+y,data=don)
pres2 <- predict(res2,dongrid)
## loading data
load("res.rda")
load("pres.rda")
test_that("modelling results with kernel default df", {
   expect_true(all(abs(res$fitted-res2$fitted)<1e-08))
   expect_true(all(abs(res$iter-res2$iter)<1e-08))
   expect_true(all(abs(res$initialdf-res2$initialdf)<1e-08))
   expect_true(all(abs(res$finaldf-res2$finaldf)<1e-08))
   expect_true(all(abs(res$bandwidth-res2$bandwidth)<1e-08))
   expect_true(all(abs(res$criteria-res2$criteria)<1e-08))
    expect_true(all(abs(pres-pres2)<1e-08))
})
## TPS
res2 <- ibr(z~x+y,data=don,smoother="tps")
pres2 <- predict(res2,dongrid)
## loading data
load("restps.rda")
load("prestps.rda")
test_that("modelling results with TPS splines default df", {
   expect_true(all(abs(res$fitted-res2$fitted)<1e-08))
   expect_true(all(abs(res$iter-res2$iter)<1e-08))
   expect_true(all(abs(res$initialdf-res2$initialdf)<1e-08))
   expect_true(all(abs(res$finaldf-res2$finaldf)<1e-08))
   expect_true(all(abs(res$bandwidth-res2$bandwidth)<1e-08))
   expect_true(all(abs(res$criteria-res2$criteria)<1e-08))
    expect_true(all(abs(pres-pres2)<1e-08))
})
## LOWRANK TPS
res2 <- ibr(z~x+y,data=don,smoother="lrtps",rank=80)
## loading data
load("reslrtps.rda")
test_that("modelling results with lowrank TPS splines default df and rank=80", {
   expect_true(all(abs(res$fitted-res2$fitted)<1e-06))
   expect_true(all(abs(res$iter-res2$iter)<1e-06))
})
