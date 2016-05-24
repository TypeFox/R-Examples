context("R Random Generator")

test_that("RRG sample works", {
  set.seed(100)
  for (i in 1:100) {
    smpl <- test_RRG_sample() #nolint
    expect_true(is.numeric(smpl) & is.finite(smpl))
    expect_that(smpl, is_less_than(1))
    expect_that(smpl, is_more_than(0))
  }
})

test_that("RRG sampleUnitExpo works", {
  set.seed(100)
  for (i in 1:100) {
    smpl <- test_RRG_sampleUnitExpo() #nolint
    expect_true(is.numeric(smpl) & is.finite(smpl))
    expect_that(smpl, is_more_than(0))
  }
})

test_that("RRG sampleExpoExpoLimit works", {
  set.seed(100)
  for (i in 1:100) {
    smpl <- test_RRG_sampleExpoExpoLimit(1, 0, 2) #nolint
    expect_true(is.numeric(smpl) & is.finite(smpl))
    expect_true(smpl > 0 | smpl == -1)
    expect_that(smpl, is_less_than(2))
  }
})
