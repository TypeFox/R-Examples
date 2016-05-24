# Test the grid optimization function


# Options
# expect_that(x, is_true())   expect_true(x)
# expect_that(x, is_false())   expect_false(x)
# expect_that(x, is_a(y)) 	expect_is(x, y)
# expect_that(x, equals(y)) 	expect_equal(x, y)
# expect_that(x, is_equivalent_to(y)) 	expect_equivalent(x, y)
# expect_that(x, is_identical_to(y)) 	expect_identical(x, y)
# expect_that(x, matches(y)) 	expect_match(x, y)
# expect_that(x, prints_text(y)) 	expect_output(x, y)
# expect_that(x, shows_message(y)) 	expect_message(x, y)
# expect_that(x, gives_warning(y)) 	expect_warning(x, y)
# expect_that(x, throws_error(y)) 	expect_error(x, y)


context( "Optimization works" )

test_that( "Optimization test 1: simple quadratic with simple inequality", {
  set.seed(100)
  f <- function(x) x^2
  b <- function(x) x>5
  res <- sequentialGridSearch( f=f, constraint=b, bounds=c(0,20), verbose=TRUE )
  expect_true( res > 4.999 & res < 5.001 )
} )

test_that( "Optimization test 2: Search for lower CI bound", {
  set.seed(100)
  lprb <- sequentialGridSearch( 
    f=identity, # We just want to minimize pr
    constraint=function(probs,...) vapply( probs, FUN=medianConsistentlyOne, FUN.VALUE=NA, ... ),
    bounds=c(0,1), 
    verbose=FALSE,
    size=100, R=50000,
    nEach=40, shrink=10, warn=FALSE
  )
  expect_true( lprb > 0.993 & lprb < 0.994 )
} )