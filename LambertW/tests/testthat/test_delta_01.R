
delta.v <- c(0, 0.1, 0.3)

test_that("delta_01 returns numeric with correct names and correct values", {
  for (dd in delta.v) {
    params.dd <- delta_01(dd)
    expect_identical(names(params.dd), 
                     c("mu_x", "sigma_x", "gamma", "delta", "alpha")) 
    expect_equivalent(params.dd["delta"], dd)
    expect_equivalent(params.dd["mu_x"], 0)
    expect_equivalent(params.dd["gamma"], 0)
    expect_equivalent(params.dd["alpha"], 1)
  }
})

test_that("delta_01 returns NA if delta is too large", {
  
  params.01.half <- delta_01(0.5)
  expect_true(is.na(params.01.half["sigma_x"]))
  expect_equivalent(params.01.half["mu_x"], 0)
  
  params.01.one <- delta_01(1)
  expect_true(all(is.na(params.01.one)[c("mu_x", "sigma_x")]))
})

test_that("delta_01 arguments of mu and sigma work", {
  params.15.half <- delta_01(0.1, mu.y = 1, sigma.y = 5)
  
  expect_true(params.15.half["sigma_x"] < 5)
  expect_equivalent(params.15.half["mu_x"], 1)
  
})

test_that("throws error for non Gaussian distribution", {
  expect_error(delta_01(0, distname = "exp"))
})
