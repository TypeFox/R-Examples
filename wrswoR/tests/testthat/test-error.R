context("error")

test_that("sample larger than population", {
  for (funcname in funcnames) {
    expect_error(funcs[[funcname]](1, 2, 1), "larger than the population", label = funcname)
  }
})

test_that("number of probabilities", {
  for (funcname in funcnames) {
    expect_error(funcs[[funcname]](3, 2, 1:2), "incorrect number of probabilities", label = funcname)
    expect_error(funcs[[funcname]](3, 2, 1:4), "incorrect number of probabilities", label = funcname)
    expect_error(funcs[[funcname]](4, 3, 1:5), "incorrect number of probabilities", label = funcname)
  }
})
