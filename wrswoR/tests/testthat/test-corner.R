context("corner")

test_that("zero-length output", {
  for (funcname in funcnames) {
    expect_equal(funcs[[funcname]](10, 0, 1:10), integer(), label = funcname)
  }
})

test_that("uniform probability", {
  n <- 13
  prob <- rep(1, n)
  for (funcname in funcnames) {
    expect_equal(sort(funcs[[funcname]](n, n, prob)), seq_len(n), label = funcname)
  }
})

test_that("heavily skewed probability", {
  prob <- c(rep(1, 10), 1e10)
  for (funcname in funcnames) {
    expect_equal(funcs[[funcname]](length(prob), 1, prob), length(prob), label = funcname)
  }
})

test_that("heavily skewed probability (reversed)", {
  prob <- c(1e10, rep(1, 10))
  for (funcname in funcnames) {
    expect_equal(funcs[[funcname]](length(prob), 1, prob), 1, label = funcname)
  }
})
