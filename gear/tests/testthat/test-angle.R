coords1 = matrix(0, nrow = 8, ncol = 2)
coords2 = cbind(c(2, 2, 0, -2, -2, -2, 0, 2), c(0, 2, 2, 2, 0, -2, -2, -2))

test_that("check accuracy of angle2d function", {
  expect_that(angle2d(coords1, coords2), equals(c(0,  45,  90, 135,   0,  45,  90, 135)))
})
