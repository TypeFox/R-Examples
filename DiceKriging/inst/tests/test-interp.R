d <- 2; n <- 16
design.fact <- expand.grid(x1=seq(0,1,length=4), x2=seq(0,1,length=4))
y <- apply(design.fact, 1, branin) 

# kriging model 1 : matern5_2 covariance structure, no trend
m1 <- km(design=design.fact, response=y, 
         coef.trend=130, coef.cov=c(0.3, 0.8), coef.var=10000)
# with nugget: should still interpolate (a difference between noisy observations) 
m1Nugget <- km(design=design.fact, response=y, nugget = 1000, 
         coef.trend=130, coef.cov=c(0.3, 0.8), coef.var=10000)


p <- predict(m1, newdata=design.fact, type="UK")
pNugget <- predict(m1Nugget, newdata=design.fact, type="UK")

precision <- 1e-10  # the following tests should work with it, since the computations are analytical
test_that(desc="Kriging mean (no nugget), on the design points", 
          expect_true(max(abs(p$mean - y)) < precision))
test_that(desc="Kriging mean with nugget, on the design points", 
          expect_true(max(abs(pNugget$mean - y)) < precision))


