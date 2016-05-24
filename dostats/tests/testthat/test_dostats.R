library(plyr)

test_that("Testing dostats",{
  expect_is(dostats(1:3, mean, sd), 'data.frame')
  expect_identical(dostats(1:3, mean, sd), data.frame(mean=2,sd=1))
  expect_identical(names(dostats(1:3, mean, sd, N=length)), c("mean", 'sd', 'N'))  
  expect_error(dostats(NA, mean))
  if(require(plyr)){
    expect_warning(ldply(iris, dostats, mean), info = "Factor passed to mean")
    expect_error(ldply(iris, dostats, median), info = "factor passed to median")
    nsd <- ldply(iris, numeric.stats, mean, median)
    expect_is(nsd, 'data.frame')
    expect_equal(dim(nsd), c(4,3), info="Checking dimentions of resulting data frame.")
    fsd <- ldply(iris, numeric.stats, mean, median)
    expect_is(fsd, 'data.frame')
    expect_equal(dim(fsd), c(4,3), info="Checking dimentions of resulting data frame.")
    rm(nsd, fsd)
  }
})
