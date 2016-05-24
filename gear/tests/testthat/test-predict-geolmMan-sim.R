set.seed(10)

# generate response
y = rnorm(10)
# generate coordinates
x1 = runif(10); x2 = runif(10)

# data frame for observed data
data = data.frame(y, x1, x2)
coords = cbind(x1, x2)
d = as.matrix(dist(coords))
psill = 2 # partial sill
r = 4 # range parameter
evar = .1 # error variance
fvar = .1 # add finescale variance
# one can't generally distinguish between evar and fvar, but
# this is done for illustration purposes

# manually specify a an expoential covariance model 
v = psill * exp(-d/r) + (evar + fvar) * diag(10)

cmod_man = cmod.man(v = v, evar = evar)

# geolm for universal kriging
gearmod_uk = geolm(y ~ x1 + x2, data = data,
                 coordnames = c("x1", "x2"),
                 cmod = cmod_man)
# newdata must have columns with prediction coordinates
# add 5 unsampled sites to sampled sites
newdata = data.frame(x1 = c(x1, runif(5)), x2 = c(x2, runif(5)))
newcoords = newdata[,c("x1", "x2")]
# create vop and vp using distances
dop = sp::spDists(as.matrix(coords), as.matrix(newcoords))
dp = as.matrix(dist(newcoords))

vop = psill * exp(-dop/r) + fvar * (dop == 0)
vp = psill * exp(-dp/r) + fvar * diag(nrow(newcoords))

# prediction for universal kriging, with conditional simulation
set.seed(11)
pred_uk_man = predict(gearmod_uk, newdata, vop = vop, vp = vp, 
                      dmethod = "svd", nsim = 10, sp = FALSE)

cmod_std = cmod.std("exponential", psill = psill, r = r, evar = evar, fvar = fvar)

set.seed(11)
gearmod_uk2 = geolm(y ~ x1 + x2, data = data,
                  coordnames = c("x1", "x2"),
                  cmod = cmod_std)
pred_uk_std = predict(gearmod_uk2, newdata, nsim = 10, dmethod = "svd", sp = FALSE)

test_that("predict.geolmMan uk calculations are correct", {
  expect_true(max(abs(range(pred_uk_man$pred - pred_uk_std$pred))) < 1e-10)
  expect_true(max(abs(range(pred_uk_man$mspe - pred_uk_std$mspe))) < 1e-10)
  expect_true(max(abs(range(pred_uk_man$sim - pred_uk_std$sim))) < 1e-7)
})

# do the same thing for simple kriging
gearmod_sk = geolm(y ~ x1 + x2, data = data,
                   coordnames = c("x1", "x2"),
                   cmod = cmod_man, mu = 2)

# prediction for simple kriging, with conditional simulation
set.seed(12)
pred_sk_man = predict(gearmod_uk, newdata, nsim = 2, vop = vop, vp = vp, 
                      dmethod = "svd", sp = FALSE)

# compare to std approach
gearmod_sk2 = geolm(y ~ x1 + x2, data = data,
                    coordnames = c("x1", "x2"),
                    cmod = cmod_std, mu = 2)
set.seed(12)
pred_sk_std = predict(gearmod_uk2, newdata, nsim = 2, dmethod = "svd", sp = FALSE)

test_that("predict.geolmMan sk calculations are correct", {
  expect_true(max(abs(range(pred_sk_man$pred - pred_sk_std$pred))) < 1e-10)
  expect_true(max(abs(range(pred_sk_man$mspe - pred_sk_std$mspe))) < 1e-10)
  expect_true(max(abs(range(pred_sk_man$sim - pred_sk_std$sim))) < 1e-7)
})

