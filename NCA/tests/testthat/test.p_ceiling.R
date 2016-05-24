test_that("str_length is number of characters", {
  scope = list(x.low=1, x.high=6, y.low=2, y.high=6)

  expect_equal(p_ceiling(scope, 0,    8), 0)
  expect_equal(p_ceiling(scope, -2/3, 8), NaN) # 3
  expect_equal(p_ceiling(scope, -4/3, 8), NaN) # 12

  expect_equal(p_ceiling(scope, 2/3,  4), 4/3)
  expect_equal(p_ceiling(scope, 0,    4), 10)
  expect_equal(p_ceiling(scope, -2/3, 4), NaN) # 4/3

  expect_equal(p_ceiling(scope, 4/3,  0), 8)
  expect_equal(p_ceiling(scope, 2/3,  0), 17)
  expect_equal(p_ceiling(scope, 0,    0), 0)
})
