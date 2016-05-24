# define context
context("leontief")

# load test data
data(leather)

# leontief decomposition
# n.b. using default method (Leontief)
# with default post-multiplication (exports)
l <- decomp(inter,
            final,
            countries,
            industries,
            out)

# test output format (i.e. structure not numbers)
test_that("output size matches", {
  expect_equal( length(l), 5 )
  expect_equal( dim(l)[1], 81 )
})

test_that("output format matches", {
  expect_match( typeof(l[,5]), "double" )
})

# test output content (i.e. numbers)
test_that("output matches", {
  expect_equal( l[1, 5],  28.52278, tolerance = .002 )
  expect_equal( l[81, 5], 34.74381, tolerance = .002 )
})


context("leontief-output")

# leontief decomposition
lo <- decomp(inter,
             final,
             countries,
             industries,
             out,
             method = "leontief",
             post = "output")

test_that("output size matches", {
  expect_equal( length(lo), 5 )
  expect_equal( dim(lo)[1], 81 )
})

test_that("output format matches", {
  expect_match(typeof( lo[,5]), "double" )
})

# test output content (i.e. numbers)
test_that("output matches", {
  expect_equal( lo[1, 5],  66.75361799, tolerance = .002 )
  expect_equal( lo[81, 5], 96.78316785, tolerance = .002 )
})
