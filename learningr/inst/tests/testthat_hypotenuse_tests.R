expect_equal(learningr::hypotenuse(3, 4), 5)
expect_error(learningr::hypotenuse())
expect_equal(
  learningr::hypotenuse(1e-300, 1e-300), 
  sqrt(2) * 1e-300, 
  tol = 1e-305
)
expect_equal(learningr::hypotenuse(1e300, 1e300), sqrt(2) * 1e300)
