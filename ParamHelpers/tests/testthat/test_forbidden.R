context("forbidden")


test_that("forbidden", {

  ps = makeParamSet(
    makeNumericParam("x", lower = 1, upper = 5),
    makeIntegerParam("y", lower = 1, upper = 10),
    makeLogicalParam("z"),
    forbidden = quote(x > 2)
  )

  expect_true(!isForbidden(ps, list(x = 1, y = 7, z = TRUE)))
  expect_true( isForbidden(ps, list(x = 3, y = 7, z = TRUE)))
  expect_true( isFeasible(ps,  list(x = 1, y = 7, z = TRUE)))
  expect_true(!isFeasible(ps,  list(x = 3, y = 7, z = TRUE)))

  xs = sampleValues(ps, 1000)
  fb = sapply(xs, isForbidden, par.set = ps)
  expect_true(!any(fb))
  ok = sapply(xs, isFeasible, par = ps)
  expect_true(all(ok))

  d = generateGridDesign(ps, resolution = 10)
  expect_true(all(d$x <= 2))

  d = generateDesign(1000, ps)
  expect_true(all(d$x <= 2))
})


