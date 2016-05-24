set.seed(10)

# generate response
y = rnorm(10)
# generate coordinates
x1 = runif(10); x2 = runif(10)

# data frame for observed data
data = data.frame(y, x1, x2)
coords = cbind(x1, x2)
psill = 2 # partial sill
r = 4 # range parameter
evar = .1 # error variance
fvar = .1 # add finescale variance
# one can't generally distinguish between evar and fvar, but
# this is done for illustration purposes

cmod_std = cmod.std("exponential", psill = psill, r = r, 
                    evar = evar, fvar = fvar)

cmod_std2 = cmod.std("exponential", psill = psill + 1, r = r + .5, 
                    evar = evar + .01, fvar = fvar)

# check geolm update for universal kriging
gear1 = geolm(y ~ x1 + x2, data = data,
                 coordnames = c("x1", "x2"),
                 cmod = cmod_std)
gear2 = geolm(y ~ x1 + x2, data = data,
              coordnames = c("x1", "x2"),
              cmod = cmod_std2)
gear2b = update(gear1, cmod_std2)

test_that("update.geolmStd uk calculations are correct", {
  expect_true(identical(gear2, gear2b))
})

# geolm for simple kriging
gear1 = geolm(y ~ x1 + x2, data = data,
              coordnames = c("x1", "x2"),
              cmod = cmod_std, mu = 2)
gear2 = geolm(y ~ x1 + x2, data = data,
              coordnames = c("x1", "x2"),
              cmod = cmod_std2, mu = 2)
gear2b = update(gear1, cmod_std2)

test_that("update.geolmStd sk calculations are correct", {
  expect_true(identical(gear2, gear2b))
})



