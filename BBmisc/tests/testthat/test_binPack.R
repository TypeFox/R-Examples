context("binPack")

test_that("binPack", {
  x = 1:10
  res = binPack(x, 11L)
  sums = sapply(split(x, res), sum)
  expect_true(is.integer(res))
  expect_equal(sums, setNames(rep(11L, 5L), 1:5))
  expect_true(all(sums <= 11L))

  x = sample(seq(from=0, to=1, by=0.01))
  res = binPack(x, 1)
  sums = sapply(split(x, res), sum)
  expect_true(is.integer(res))
  expect_true(all(head(sums, 50) == 1))

  x = runif(20)
  res = binPack(x, 1)
  sums = sapply(split(x, res), sum)
  expect_true(is.integer(res))
  expect_true(all(sums < 1))

  x = runif(5)
  res = binPack(x, Inf)
  expect_true(is.integer(res))
  expect_true(length(res) == 5 && all(as.numeric(res) == 1))

  expect_error(binPack(c(-5, 3)))
  expect_error(binPack(c(1, 100), 10))
  expect_error(binPack(c(1, Inf), 1))
})
