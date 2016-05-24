n = 10
coords = matrix(runif(n * 2), nrow = n)
d = as.matrix(dist(coords))

# check geoR
mod = c("exponential", "matern", "spherical", "gaussian")
evar = c(0, .5, 0, .5)
fvar = c(0, 0, .5, .5)
r = runif(4)
par3 = runif(4, 0, 2)
psill = rgamma(4, 2, 2)

compare = numeric(length(mod)*4)
count = 0
for(i in 1:length(mod))
{
  for(j in 1:length(evar))
  {
    count = count + 1
    A = varcov.spatial(coords, cov.model = mod[i], 
                       kappa = par3[j],
                       nugget = (evar[j] + fvar[j]), 
                       cov.pars = c(psill[j], r[j]))
    cmod = cmod.std(model = mod[i], 
                   par3 = par3[j],
                   psill = psill[j], 
                   r = r[j], 
                   evar = evar[j],
                   fvar = fvar[j])
    B = eval.cmod(cmod, d)
    compare[count] = max(abs(range(A$varcov - B)))
  }
}

test_that("eval.cmod.cmodStd is accurate (geoR)", {
  expect_true(max(compare) < 1e-10)
})

# comparison with spam

mod = c("exponential", "spherical", "wendland1", "wendland2", "wu1", "wu2", "wu3")
evar = c(0, .5, 0, .5)
fvar = c(0, 0, .5, .5)
r = runif(4)
par3 = runif(4, 0, 2)
psill = rgamma(4, 2, 2)

compare = numeric(length(mod)*4)
count = 0
for(i in 1:length(mod))
{
  for(j in 1:length(evar))
  {
    count = count + 1
    A = covmat(d, theta = c(r[j], psill[j], (evar[j]+fvar[j])),
               type = mod[i])
    cmod = cmod.std(model = mod[i], 
                   par3 = par3[j],
                   psill = psill[j], 
                   r = r[j], 
                   evar = evar[j],
                   fvar = fvar[j])
    B = eval.cmod(cmod, d)
    compare[count] = max(abs(range(A - B)))
  }
}

test_that("eval.cmod.cmodStd is accurate (spam)", {
  expect_true(max(compare) < 1e-10)
})

# comparison with spam (matern)
mod = c("matern")
evar = c(0, .5, 0, .5)
fvar = c(0, 0, .5, .5)
r = runif(4)
par3 = runif(4, 0, 2)
psill = rgamma(4, 2, 2)

compare = numeric(length(mod)*4)
count = 0
for(i in 1:length(mod))
{
  for(j in 1:length(evar))
  {
    count = count + 1
    A = covmat(d, theta = c(r[j], psill[j], par3[j], (evar[j]+fvar[j])),
               type = mod[i])
    cmod = cmod.std(model = mod[i], 
                   par3 = par3[j],
                   psill = psill[j], 
                   r = r[j], 
                   evar = evar[j],
                   fvar = fvar[j])
    B = eval.cmod(cmod, d)
    compare[count] = max(abs(range(A - B)))
  }
}

test_that("eval.cmod.cmodStd matern is accurate (spam)", {
  expect_true(max(compare) < 1e-10)
})



