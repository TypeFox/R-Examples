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
cmod_std0 = cmod_std
cmod_std0$evar = 0
v1 = eval.cmod(cmod_std0, 
                     d = sp::spDists(coords)) +
           evar * diag(10)
cmod_man = cmod.man(v = v1, evar = cmod_std$evar)
cmod_std2 = cmod.std("exponential", psill = psill + 1, r = r + .5, 
                    evar = evar + .01, fvar = fvar)
cmod_std20 = cmod_std
cmod_std20$evar = 0
v2 = eval.cmod(cmod_std20, 
                     d = sp::spDists(coords)) +
              cmod_std2$evar * diag(10)
cmod_man2 = cmod.man(v = v2, evar = cmod_std2$evar)

# check geolm update for universal kriging
gear1 = geolm(y ~ x1 + x2, data = data,
                 coordnames = c("x1", "x2"),
                 cmod = cmod_man)
gear2 = geolm(y ~ x1 + x2, data = data,
              coordnames = c("x1", "x2"),
              cmod = cmod_man2)
gear2b = update(gear1, cmod_man2)

test_that("update.geolmStd uk calculations are correct", {
  expect_true(identical(gear2, gear2b))
})

# geolm for simple kriging
gear1 = geolm(y ~ x1 + x2, data = data,
              coordnames = c("x1", "x2"),
              cmod = cmod_man, mu = 2)
gear2 = geolm(y ~ x1 + x2, data = data,
              coordnames = c("x1", "x2"),
              cmod = cmod_man2, mu = 2)
gear2b = update(gear1, cmod = cmod_man2)

test_that("update.geolmStd sk calculations are correct", {
  expect_true(identical(gear2, gear2b))
})



