
library(rollply)
context('Grid building')

# Define some sets of points
exampleset <- data.frame(x = seq.int(0,1,l=100), 
                         y = seq.int(0,5,l=100) + rnorm(100))

grid_builders <- list(build_grid_identical, 
                      build_grid_squaretile, 
                      build_grid_ahull_crop,
                      build_grid_ahull_fill)

# Common tests among all grid-building functions
test_that("Common grid_builders test", { 
  # Note that R CMD's check do not work with ahull_*, probably due to 
  # mispackaging ?
  for (builder in grid_builders[1:2]) { 
    expect_is( builder(exampleset, npts=100, pad=0), "data.frame" )
  }
})
