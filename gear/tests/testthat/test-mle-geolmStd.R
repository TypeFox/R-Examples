set.seed(1)
library(geoR)
fields <- grf(100, cov.pars = c(1, 1), nugget = 0.1, cov.model = "exponential", 
              kappa = 0.5, messages = FALSE)
geoR_fit = likfit(fields, trend = "1st", ini.cov.pars = c(1, 1), 
                  nugget = 0.1, lik.method = "ML")

cmod1 = cmod.std("exponential", psill = 1, r = 1, evar = 0.1)
data = data.frame(y = fields$data, x1 = fields$coords[,1], x2 = fields$coords[,2])
object = geolm(y ~ x1 + x2, data = data, coordnames = c("x1", "x2"), 
               cmod = cmod1)
gear_fit = mle(object, reml = FALSE)

# make sure the results are similar
test_that("mle.geolmStd calculations are correct", {
  expect_true(abs(gear_fit$cmod_evar0$psill - geoR_fit$cov.pars[1]) < 1e-3)
  expect_true(abs(gear_fit$cmod_evar0$r - geoR_fit$cov.pars[2]) < 1e-3)
  expect_true(abs(gear_fit$evar - geoR_fit$nugget) < 1e-3)
  expect_true(abs(gear_fit$loglik - geoR_fit$loglik) < 1e-3)
})

geoR_fit_reml = likfit(fields, trend = "1st", ini.cov.pars = c(1.6, 3.7), 
                  nugget = 0.1, lik.method = "REML")
gear_fit_reml = mle(object, reml = TRUE)

test_that("mle.geolmStd calculations are correct", {
   expect_true(abs(gear_fit_reml$cmod_evar0$psill - geoR_fit_reml$cov.pars[1]) < 2e-2)
   expect_true(abs(gear_fit_reml$cmod_evar0$r - geoR_fit_reml$cov.pars[2]) < 4e-2)
   expect_true(abs(gear_fit_reml$evar - geoR_fit_reml$nugget) < 1e-3)
#   expect_true(abs(gear_fit_reml$loglik - geoR_fit_reml$loglik) < 1e-3)
})

# compare to reml likelihood manually (p. 117 of Model Based Geostatistics)
n = length(object$y)
# manually create covariance matrix
v = gear_fit_reml$cmod_evar0$psill * exp(-as.matrix(dist(object$coords))/
                                           gear_fit_reml$cmod_evar0$r) + 
  gear_fit_reml$evar * diag(n)
x = object$x
yhat = x %*% gear_fit_reml$coeff
resid = object$y - yhat
minus2_reml_ll = (n * log(2 * pi) + 
                   determinant(gear_fit_reml$v, logarithm = TRUE)$mod + 
          determinant(crossprod(x, solve(v, x)), logarithm = TRUE)$mod +
          crossprod(resid, solve(v, resid))[1, 1])

# make sure the results are similar
test_that("mle.geolmStd reml objective function correct", {
  expect_true(abs(gear_fit_reml$optimx$value - minus2_reml_ll) < 1e-3)
})

# par = c(1, .1)
# weights = rep(1, 100)
# x = cbind(1, fields$coords)
# y = fields$data + x %*% c(1, 2, 1)
# scmod = cmod.std("exponential", psill = 1, r = 1, evar = 0)
# d = as.matrix(dist(fields$coords))
# cmod = cmod.std("exponential", psill = 1, r = 1, evar = 0.1)
# data = data.frame(y = y, x1 = fields$coords[,1], x2 = fields$coords[,2])
# object = geolm(y ~ x1 + x2, data = data, coordnames = c("x1", "x2"), 
#                cmod = cmod)
# 
# nugget = "e"
# reml = FALSE


