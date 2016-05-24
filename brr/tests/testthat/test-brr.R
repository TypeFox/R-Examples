context("Brr")

test_that("Test prior function", {
  expect_equal(prior(list(a=2,b=3)), "semi-informative")
  expect_equal(prior(list(a=2,b=3,c=0.5,d=0)), "semi-informative")
  expect_equal(prior(list(a=2,b=3,c=NULL,d=NULL)), "semi-informative")
  expect_equal(prior(list()), "non-informative")
  expect_equal(prior(list(c=NULL)), "non-informative")  
  expect_equal(prior(list(a=2,b=3,c=1,d=2)), "informative")
  expect_error(prior(list(a=2)), "Invalid combination of parameters a,b,c,d.")
  expect_error(prior(list(a=2, b=0)), "Invalid combination of parameters a,b,c,d.")
  expect_error(prior(list(a=2, b=1, c=2)), "Invalid combination of parameters a,b,c,d.")
  expect_error(prior(list(a=2, b=1, c=0)), "Invalid combination of parameters a,b,c,d.")
  expect_error(prior(list(a=2, b=1, c=0.5)), "Invalid combination of parameters a,b,c,d.")
  expect_error(prior(list(a=2, b=1, c=0, d=4)), "Invalid combination of parameters a,b,c,d.")
  expect_error(prior(list(a=NULL, b=0)), "Invalid combination of parameters a,b,c,d.")
  expect_error(prior(list(c=2, d=3)), "Invalid combination of parameters a,b,c,d.")
  expect_error(prior(list(c=0.5, d=NULL)), "Invalid combination of parameters a,b,c,d.")
})

test_that("Test sprior function", {
  expect_is(sprior(Brr(a=2,b=3), "mu"), "list")
  expect_identical(sprior(Brr(a=2,b=3), "mu"), sprior_mu(a=2,b=3))
})

test_that("Test Brr function", {
  model <- Brr(a=2, b=4)
  expect_that(
    identical(model(), list(a=2,b=4)), 
    is_true()
  )
  expect_error(model(hello=10))
  model <- model(a=3)
  expect_that(
    identical(model(), list(a=3,b=4)), 
    is_true()
  )
  model <- model(c=4)
  expect_that(
    identical(model()[c("a","b","c")], list(a=3,b=4,c=4)), 
    is_true()
  )
})

test_that("Test dprior function", {
  model <- Brr(a=2, b=4)
  expect_that(
    all(dprior(model, "mu", 1:3) == dprior_mu(mu=1:3, a=2, b=4)), 
    is_true()
  )
  model <- Brr(a=2, b=4, c=3, d=10, T=10)
  expect_that(
    all(dprior(model, "x", 1:3) == dprior_x(x=1:3, a=2, b=4, c=3, d=10, T=10)), 
    is_true()
  )
})

test_that("Test dpost function", {
  # mu
  model <- Brr(a=2, b=4)
  expect_error(
    dpost(model, "mu", 1)
  )
  model <- model(c=3, d=4, x=5, y=6, T=8)
  expect_that(
    # car brr:::prior(model())==semi-info car il manque S 
    all(dpost(model, "mu", 1:3) == dpost_mu(mu=1:3, a=2, b=4, c=3, d=4, x=5, y=6, T=8)), 
    is_true()
  )  
  ## the same :
  expect_equal(dpost(model, "mu", 1:3),  dpost_mu(mu=1:3, a=2, b=4, c=3, d=4, x=5, y=6, T=8))
  # phi
  expect_error(
    dpost(model, "phi", 1)
  )
  model <- model(S=10)
  expect_that(
    all(dpost(model, "phi", 1:3) == dpost_phi(phi=1:3, a=2, b=4, c=3, d=4, x=5, y=6, S=10, T=8)), 
    is_true()
  )  
  # x 
  expect_error(dpost(model, "x", 1:3), "Missing parameters. You must supply Snew, a, c, d, x, y, S.")
  model <- model(Snew=12)
  expect_that(
    all(dpost(model, "x", 1:3) == dpost_x(xnew=1:3, a=2, c=3, d=4, S=10, x=5, y=6, Snew=12)), 
    is_true()
  )
})

test_that("Test plot.brr",{
  model <- Brr(a=2, b=4, c=3, d=10, S=11, T=10, x=2, y=3)
  expect_equal(plot(model, dpost(lambda)), invisible())
  expect_equal(plot(model, dprior(lambda)), invisible())
  expect_error(plot(model, ppost(lambda)), "ppost_lambda does not exist in brr package")
  expect_error(plot(model, pprior(lambda)))
  model <- model(c=2.9)
  expect_equal(plot(model, pprior(lambda)), invisible())
})